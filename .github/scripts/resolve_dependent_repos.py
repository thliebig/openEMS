#!/usr/bin/env python3
#
# Copyright (C) 2026 Yifeng Li <tomli@tomli.me>
#
# Permission to use, copy, modify, and/or distribute this software for
# any purpose with or without fee is hereby granted.
#
# THE SOFTWARE IS PROVIDED “AS IS” AND THE AUTHOR DISCLAIMS ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE
# FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY
# DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN
# AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
# OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
#
# To test it locally, you can create a shell script:
#
# #!/bin/bash
# export GITHUB_API_URL="https://api.github.com"
# export GITHUB_REPOSITORY_OWNER="thliebig"
# export GITHUB_EVENT_NAME="pull_request"
# export GITHUB_OUTPUT="/dev/null"
# export GITHUB_FORK_OWNER="find the username who has an open pull request"
# export GITHUB_HEAD_REF="the source branch of the same pull request"
# export GITHUB_TOKEN="login is optional, but can bypass low API rate limit"
#
# Note that all variables but GITHUB_TOKEN and GITHUB_FORK_OWNER are default
# variables on GitHub Actions, while GITHUB_TOKEN and GITHUB_FORK_OWNER is are
# non-standard variables set manually in the ".yml" file.
#
# or --project openEMS, --project QCSXCAD, --project AppCSXCAD
# python3 resolve_dependent_repos.py --project CSXCAD

import argparse
import os
from time import sleep

import json
import urllib.request


PROJECT_FOUNDER = "thliebig"


def print_log(loglevel, string, title=None):
    if loglevel not in ["debug", "notice", "warning", "error"]:
        raise ValueError("Unsupported GitHub Actions loglevel: %s" % loglevel)

    parameters = ""
    if title:
        parameters = " title=%s" % title

    print("::%s%s::%s" % (loglevel, parameters, string))


def getenv(var):
    value = os.getenv(var)
    if not value:
        errormsg = "Unable to determine GitHub Actions variable $%s!" % var
        print_log("error", errormsg)
        raise ValueError(errormsg)
    else:
        return value


def get_default_repo_dependency(project):
    REPO = {
        "fparser": {
            "owner": PROJECT_FOUNDER,
            "name": "fparser",
            "branch": "master",
        },
        "CSXCAD": {
            "owner": PROJECT_FOUNDER,
            "name": "CSXCAD",
            "branch": "master"
        },
        "openEMS": {
            "owner": PROJECT_FOUNDER,
            "name": "openEMS",
            "branch": "master"
        },
        "QCSXCAD": {
            "owner": PROJECT_FOUNDER,
            "name": "QCSXCAD",
            "branch": "master"
        },
        "AppCSXCAD": {
            "owner": PROJECT_FOUNDER,
            "name": "AppCSXCAD",
            "branch": "master"
        },
    }

    # We can check all repos, but it wastes GitHub API calls.
    # Especially for local testing without logging it, the quota
    # is only 60 calls/hr. You can burn it out just after a few
    # try. It also spams useless NOTICE messages per repo to the
    # GitHub Actions Summary.
    #
    # So we only check repos actually used by "project".
    PER_PROJECT_DEPENDENCY = {
        "CSXCAD": [REPO["fparser"]],
        "openEMS": [REPO["fparser"], REPO["CSXCAD"]],
        "QCSXCAD": [REPO["fparser"], REPO["CSXCAD"]],
        "AppCSXCAD": [REPO["fparser"], REPO["CSXCAD"], REPO["QCSXCAD"]],
    }
    return PER_PROJECT_DEPENDENCY[project]


def determine_repo_dependency(project):
    github_event_name = getenv("GITHUB_EVENT_NAME")
    if github_event_name not in ["pull_request", "push"]:
        print_log("error", "Unable to determine $GITHUB_EVENT_NAME!")
        return get_default_repo_dependency(project)

    if github_event_name == "pull_request":
        src_branch = getenv("GITHUB_HEAD_REF")
    elif github_event_name == "push":
        src_branch = getenv("GITHUB_REF_NAME")
    else:
        assert False

    repo_dependency = get_default_repo_dependency(project)

    # If this project is a fork, we need to check whether the
    # user or organization has also forked other repos. If they
    # did, rewrite the project owner's name to the fork owner's
    # name. This allows the fork itself to become a new "upstream"
    # with correct CI/CD, since forks can also accept commits
    # and Pull Requests.
    owner = getenv("GITHUB_REPOSITORY_OWNER")
    if owner != PROJECT_FOUNDER:
        for repo in repo_dependency:
            if reponame_exists(owner, repo["name"]):
                repo["owner"] = owner

    # Furthermore, if the dependency repo has a branch with the
    # same name as the current repo branch, use that branch instead.
    # This allowing testing multi-repo features.
    for repo in repo_dependency:
        if reponame_has_branch(owner, repo["name"], src_branch):
            repo["branch"] = src_branch

            # only warn dependencies, not current repo
            current_repo = repo["name"] == getenv("GITHUB_REPOSITORY").split("/")[1]

            if not current_repo:
                print_log(
                    "warning",
                    'repo %s/%s branch %s is used instead of the default branch.' %
                    (repo["owner"], repo["name"], src_branch),
                    title='%s/%s: different branch "%s" used' %
                    (repo["owner"], repo["name"], src_branch)
                )

    # If the event is a pull request...
    if github_event_name == "pull_request":
        contributor = getenv("GITHUB_FORK_OWNER")

        # check whether the contributor has forked other dependent repos.
        for repo in repo_dependency:
            # If they didn't, use the upstream dependency.
            if not reponame_exists(contributor, repo["name"]):
                continue

            # If they did, check whether they've submitted an open pull
            # request against the upstream, and the source branch uses
            # the same name
            pr = (
                pr_by_user_with_branch(
                    repo["owner"], repo["name"], contributor, src_branch
                )
            )
            if not pr:
                continue

            # If they also did, use that Pull Request's merge commit as
            # the dependency repo.
            repo["branch"] = pr["merge_commit_sha"]
            repo["pr_url"] = pr["html_url"]
            print_log(
                "warning",
                '%s/%s: %s must be merged first before merging this commit (%s)' %
                (repo["owner"], repo["name"], pr["html_url"], pr["title"]),
                title="Must-Merge Pull Request: %s" % pr["html_url"]
            )

    return repo_dependency


def url_open(url, expect_code=None):
    DEFAULT_HEADERS = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    RETRIES = 3
    error = None

    headers = DEFAULT_HEADERS.copy()
    try:
        token = getenv("GITHUB_TOKEN")
        headers["Authorization"] = "Bearer %s" % token
    except ValueError:
        # for local testing only
        print_log(
            "warning",
            "$GITHUB_TOKEN not found! API requests are rate-limited!"
        )

    request = urllib.request.Request(url, headers=headers)

    for i in range(RETRIES):
        try:
            response = urllib.request.urlopen(request)
            return response, response.getcode()
        except urllib.error.HTTPError as e:
            if e.code == expect_code:
                return None, e.code
            else:
                error = e
        except urllib.error.URLError as e:
            error = e
        sleep(1)

    raise error


def url_check_existence(url):
    response, code = url_open(url, expect_code=404)
    if code == 404:
        return False
    elif code == 200:
        return True
    else:
        raise RuntimeError("Unexpected HTTP status code: %d" % code)


def reponame_exists(owner, repo):
    API = getenv("GITHUB_API_URL")
    url = "%s/repos/%s/%s" % (API, owner, repo)
    return url_check_existence(url)


def reponame_has_branch(owner, repo, branch):
    API = getenv("GITHUB_API_URL")
    url = "%s/repos/%s/%s/branches/%s" % (API, owner, repo, branch)
    return url_check_existence(url)


def pr_by_user_with_branch(owner, repo, user, branch):
    API = getenv("GITHUB_API_URL")
    url = "%s/repos/%s/%s/pulls?head=%s:%s" % (API, owner, repo, user, branch)
    raw_response, code = url_open(url)
    pr = json.loads(raw_response.read().decode("UTF-8"))

    if pr:
        return pr[0]
    else:
        return None


def output_repo_dependency(repo_dependency):
    for repo in repo_dependency:
        repo_value = "%s/%s" % (repo["owner"], repo["name"])
        branch_value = repo["branch"]
        if "pr_url" in repo:
            pr_value = repo["pr_url"]
        else:
            pr_value = "None"

        print_log(
            "notice", "%s, pr: %s, branch: %s" % (repo_value, pr_value, branch_value),
            title="%s repo" % repo["name"]
        )

        with open(getenv("GITHUB_OUTPUT"), "a") as output:
            output.write("%s_repo=%s\n" % (repo["name"], repo_value))
            output.write("%s_branch=%s\n" % (repo["name"], branch_value))
            output.write("%s_pr_url=%s\n" % (repo["name"], pr_value))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Resolve dependent repos')
    parser.add_argument("--project", type=str, required=True)
    args = parser.parse_args()

    repo_dependency = determine_repo_dependency(args.project)
    output_repo_dependency(repo_dependency)
