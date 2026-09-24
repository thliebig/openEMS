#!/usr/bin/env python3
"""Build price_performance.md from gpu_prices.json and gpu_performance.json.

usage: gen_price_table.py [--vast]

--vast also asks the Vast.ai CLI for the cheapest offer of each GPU, which gives
the hours of renting that a purchase has to beat.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HORN_CELLS = 2416581


def main():
    with open(os.path.join(HERE, 'gpu_performance.json')) as fh:
        perf = json.load(fh)['machines']
    with open(os.path.join(HERE, 'gpu_prices.json')) as fh:
        prices = json.load(fh)

    want_vast = '--vast' in sys.argv
    if want_vast:
        sys.path.insert(0, HERE)
        from estimate import vast_price

    rows = []
    for name, p in prices['prices'].items():
        if name not in perf:
            continue
        m = perf[name]
        horn = m['measured']['horn_fd_mcells_s']
        free = m['measured']['free_space_pml_mcells_s']
        usd = p.get('median_usd')
        # where a card still sells new for less than the used listings ask, as the
        # RTX 5090 does in this shortage, the real cost is the retail price
        new = p.get('new_usd')
        if usd and new:
            usd = min(usd, new)
        elif new and not usd:
            usd = new
        row = {
            'name': name, 'usd': usd, 'condition': p.get('condition', 'used'),
            'confidence': p.get('confidence', '?'),
            'horn': horn, 'free': free,
            'per_100': round(horn / usd * 100, 0) if usd else None,
            'note': p.get('note', ''),
        }
        if want_vast:
            row['rent'] = vast_price(name)
            row['break_even_h'] = (round(usd / row['rent']) if usd and row['rent']
                                   else None)
        rows.append(row)

    rows.sort(key=lambda r: (r['per_100'] is None, -(r['per_100'] or 0)))

    out = ['# Price and performance', '',
           'What each GPU costs against what it does for openEMS. The speeds are the',
           'measurements of [benchmarks.md](../benchmarks.md); the prices are what the',
           'cards were going for in September 2026, which is a worse kind of number',
           '(see "The prices" below).', '',
           '"Horn" is the speed on a mesh of 2.4 million cells, the size of a real job;',
           '"free space" is the speed on 27 million cells, where a big GPU stretches its',
           'legs. "Hours to break even" is the purchase price divided by the cheapest',
           'Vast.ai rental of the same GPU: rent for fewer hours than that and renting wins.',
           'A dash there means nobody was renting that GPU out when this was generated,',
           'which is itself worth knowing: the rare ones are the RTX 4070 Ti and the 3080 Ti.', '',
           'AI disclosure: collected and written up with Claude Opus 5 (Claude Code).', '',
           '## Table', '']
    cols = ['GPU', 'Price', 'Horn MCells/s', 'Free space MCells/s', 'MCells/s per $100']
    if want_vast:
        cols += ['Vast $/h', 'Hours to break even']
    out += ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    for r in rows:
        cells = [r['name'],
                 ('$%d' % r['usd']) if r['usd'] else 'unknown',
                 '%.0f' % r['horn'], '%.0f' % r['free'],
                 ('**%.0f**' % r['per_100']) if r['per_100'] else '-']
        if want_vast:
            cells += ['%.3f' % r['rent'] if r.get('rent') else '-',
                      '%d h' % r['break_even_h'] if r.get('break_even_h') else '-']
        out.append('| ' + ' | '.join(cells) + ' |')

    priced = [r for r in rows if r['usd']]
    best = priced[:3]
    def fastest_under(limit):
        under = [r for r in priced if r['usd'] <= limit]
        return max(under, key=lambda r: r['horn']) if under else None
    evens = sorted(r['break_even_h'] for r in priced if r.get('break_even_h'))

    out += ['', '## What it says', '',
            '**Value goes to the old cards.** The top of the table is %s, at %.0f, %.0f and '
            '%.0f MCells/s per $100. Every one of them is a Turing or Ampere card that the '
            'gaming market has finished with, and they run the same HIP kernels as the new '
            'ones: the GPU engine needs Pascal or later and nothing else.'
            % (', '.join(r['name'] for r in best), best[0]['per_100'], best[1]['per_100'],
               best[2]['per_100']), '']
    f300, f600 = fastest_under(300), fastest_under(600)
    if f300 and f600:
        worst = priced[-2:]   # the two worst of the value ranking, both halo cards
        out += ['**The fastest card for a budget.** Under $300 that is the %s (%.0f MCells/s '
                'on the horn); under $600 the %s (%.0f). At the other end the %s at $%d and '
                'the %s at $%d buy their speed at %.0f and %.0f MCells/s per $100, a third of '
                'what the cheap cards give.'
                % (f300['name'], f300['horn'], f600['name'], f600['horn'],
                   worst[0]['name'], worst[0]['usd'], worst[1]['name'], worst[1]['usd'],
                   worst[0]['per_100'], worst[1]['per_100']), '']
    if evens:
        out += ['**Renting usually wins.** The break even against Vast.ai runs from %d to %d '
                'hours, %d in the middle: that is months of continuous simulation before a '
                'purchase pays for itself, and the rented machine is somebody else\'s problem '
                'when it breaks. Buying makes sense for a workstation that simulates every day, '
                'or where the data cannot leave the building.'
                % (evens[0], evens[-1], evens[len(evens) // 2]), '']
    out += ['**The memory of the AI market sets these prices, not the speed of the cards.** '
            'The 24 GB cards (RTX 3090, 3090 Ti, 4090) and the 16 GB RTX 4060 Ti carry a '
            'premium for their memory that this workload does not use: openEMS needed under '
            '850 MiB of GPU memory for every benchmark here, so a cheap 8 GB card runs the '
            'same job as a 24 GB one. Buy for bandwidth, not for capacity, until a mesh no '
            'longer fits.', '']
    out += ['## The prices', '',
            prices.get('_method', ''), '']
    for name, p in sorted(prices['prices'].items()):
        if not p.get('median_usd'):
            out.append('- **%s**: no price could be verified%s'
                       % (name, (' (%s)' % p['note']) if p.get('note') else ''))
            continue
        out.append('- **%s**: $%d %s, %s, %s confidence%s'
                   % (name, p['median_usd'], p.get('condition', 'used'),
                      p.get('range_usd', 'no range'), p.get('confidence', '?'),
                      (' — ' + p['note']) if p.get('note') else ''))

    dst = os.path.join(HERE, 'price_performance.md')
    with open(dst, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('wrote %s, %d GPUs (%d priced)'
          % (dst, len(rows), sum(1 for r in rows if r['usd'])))


if __name__ == '__main__':
    main()
