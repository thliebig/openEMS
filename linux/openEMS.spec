#
# spec file for package [spectemplate]
#
# Copyright (c) 2010 SUSE LINUX Products GmbH, Nuernberg, Germany.
#
# All modifications and additions to the file contributed by third parties
# remain the property of their copyright owners, unless otherwise agreed
# upon. The license for this file, and modifications and additions to the
# file, is the same license as for the pristine package itself (unless the
# license for the pristine package is not an Open Source License, in which
# case the license is the MIT License). An "Open Source License" is a
# license that conforms to the Open Source Definition (Version 1.9)
# published by the Open Source Initiative.

# Please submit bugfixes or comments via http://bugs.opensuse.org/
#

# norootforbuild

Name:           openEMS
Version:        0.0.28
Release:        3
Summary:        Free and Open Electromagnetic Field Solver
Group:          Productivity/Scientific/Physics
License:        GPLv3
URL:            http://www.openems.de
Source0:        %{name}-%{version}.tar.bz2
Patch0:         invoke_openEMS.m.patch
Patch1:         README.patch
Patch2:         CalcNF2FF.m.patch
Patch3:         fedora17.diff
BuildRoot:      %_tmppath/%name-%version-build

# libqt4-devel is needed only to provide qmake (the Qt-libraries are not used)
# libfparser4-devel contains a static library => no runtime requirement
BuildRequires:  libqt4-devel gcc-c++ libfparser4-devel hdf5-devel tinyxml-devel CSXCAD-devel openmpi-devel vtk-devel boost-devel
Requires:       CSXCAD

# determine qt4 qmake executable
%if 0%{?fedora}
    %global qmake qmake-qt4
%else
    %global qmake qmake
%endif



%description
OpenEMS is a free and open-source electromagnetic field solver using the (EC-)FDTD method.


%prep
%setup -q
%patch0 -p1
#%if 0%{?fedora}
#%patch1 -p1
#%endif
%patch1 -p1
%patch2 -p1
%if 0%{?fedora} >= 17
%patch3 -p1
%endif

%build
ADDFLAGS="-msse" # enable at least the SSE command set (no SSE makes no sense -- way too slow)
%qmake QMAKE_CFLAGS="%optflags $ADDFLAGS" QMAKE_CXXFLAGS="%optflags $ADDFLAGS" LIB_SUFFIX="$(echo %_lib | cut -b4-)" CONFIG+=packaging openEMS.pro
make %{?_smp_mflags}
cd nf2ff
%qmake QMAKE_CFLAGS="%optflags $ADDFLAGS" QMAKE_CXXFLAGS="%optflags $ADDFLAGS" LIB_SUFFIX="$(echo %_lib | cut -b4-)" CONFIG+=packaging nf2ff.pro
make %{?_smp_mflags}
cd ..

%install
make INSTALL_ROOT=%{buildroot} install
cd nf2ff
make INSTALL_ROOT=%{buildroot} install
cd ..
find %{buildroot} -name '*.la' -exec rm -f {} ';'


%clean
rm -rf %{buildroot}


%files
%defattr(-,root,root)
%doc README COPYING TODO known_bugs known_problems
/usr/share/%{name}
/usr/bin/openEMS*
/usr/bin/nf2ff


%changelog
* Sat Jun 23 2012 Sebastian Held <sebastian.held@gmx.de> - 0.0.28-3
- display correct version
* Mon Jun 18 2012 Sebastian Held <sebastian.held@gmx.de> - 0.0.28-2
- Fedora 17 build fixes
* Sun Jun 17 2012 Sebastian Held <sebastian.held@gmx.de> - 0.0.28-1
- new upstream version
* Thu Mar 1 2012 Sebastian Held <sebastian.held@gmx.de> - 0.0.27-1
- new upstream version
* Sat Jan 21 2012 Sebastian Held <sebastian.held@gmx.de> - 0.0.26-1
- new upstream version
* Mon Jan 9 2012 Sebastian Held <sebastian.held@gmx.de> - 0.0.25-4
- added runtime dep on CSXCAD
* Thu Dec 29 2011 Sebastian Held <sebastian.held@gmx.de> - 0.0.25-3
- added README
* Sun Dec 25 2011 Sebastian Held <sebastian.held@gmx.de> - 0.0.25-2
- Fedora 16 build fix and upstream fixes
* Sun Dec 4 2011 Sebastian Held <sebastian.held@gmx.de> - 0.0.25-1
- new upstream version
* Mon Oct 3 2011 Sebastian Held <sebastian.held@gmx.de> - 0.0.24-1
- initial version
