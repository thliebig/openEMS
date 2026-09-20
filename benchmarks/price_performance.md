# Price and performance

What each GPU costs against what it does for openEMS. The speeds are the
measurements of [benchmarks.md](../benchmarks.md); the prices are what the
cards were going for in September 2026, which is a worse kind of number
(see "The prices" below).

"Horn" is the speed on a mesh of 2.4 million cells, the size of a real job;
"free space" is the speed on 27 million cells, where a big GPU stretches its
legs. "Hours to break even" is the purchase price divided by the cheapest
Vast.ai rental of the same GPU: rent for fewer hours than that and renting wins.
A dash there means nobody was renting that GPU out when this was generated,
which is itself worth knowing: the rare ones are the RTX 4070 Ti and the 3080 Ti.

AI disclosure: collected and written up with Claude Opus 5 (Claude Code).

## Table

| GPU | Price | Horn MCells/s | Free space MCells/s | MCells/s per $100 | Vast $/h | Hours to break even |
|---|---|---|---|---|---|---|
| RTX 2060 Super | $168 | 3162 | 3913 | **1882** | 0.065 | 2569 h |
| RTX 2070 | $160 | 2998 | 3772 | **1874** | 0.061 | 2606 h |
| RTX 2070 Super | $174 | 3113 | 3948 | **1789** | 0.061 | 2834 h |
| GTX 1660 Ti | $103 | 1838 | 2268 | **1784** | 0.068 | 1513 h |
| RTX 2060 | $138 | 2219 | 2674 | **1608** | 0.068 | 2028 h |
| RTX 2080 Ti | $271 | 4342 | 5803 | **1602** | 0.080 | 3383 h |
| RTX 3080 10 GB | $358 | 5416 | 7458 | **1513** | 0.121 | 2967 h |
| RTX 3070 | $260 | 3743 | 4941 | **1440** | 0.081 | 3213 h |
| RTX 3070 Ti | $309 | 4422 | 5858 | **1431** | 0.122 | 2526 h |
| RTX 3060 Ti | $248 | 3353 | 4397 | **1352** | 0.114 | 2168 h |
| GTX 1080 Ti | $155 | 1879 | 2619 | **1212** | 0.061 | 2521 h |
| RTX 3080 Ti | $487 | 5675 | 7486 | **1165** | 0.122 | 3982 h |
| RTX 5060 | $370 | 3728 | 5050 | **1008** | - | - |
| RTX 3050 | $175 | 1737 | 2205 | **993** | - | - |
| RTX 3060 | $295 | 2751 | 3649 | **933** | 0.052 | 5718 h |
| RTX 4060 | $289 | 2629 | 3536 | **910** | 0.068 | 4247 h |
| RTX 4070 | $515 | 4414 | 5894 | **857** | 0.195 | 2643 h |
| RTX 5070 | $700 | 5898 | 8160 | **843** | 0.201 | 3476 h |
| RTX 4070 Super | $550 | 4570 | 6119 | **831** | 0.200 | 2750 h |
| RTX 4070 Ti Super | $801 | 5880 | 7970 | **734** | 0.189 | 4249 h |
| RTX 4070 Ti | $688 | 4864 | 6675 | **707** | - | - |
| RTX 4080 | $905 | 6388 | 8898 | **706** | 0.268 | 3376 h |
| RTX 5070 Ti | $1050 | 7369 | 10346 | **702** | 0.255 | 4122 h |
| RTX 5060 Ti 16 GB | $640 | 4048 | 5516 | **632** | 0.116 | 5512 h |
| RTX 3090 | $1032 | 6498 | 9439 | **630** | 0.149 | 6916 h |
| RTX 4080 Super | $1095 | 6563 | 9142 | **599** | 0.219 | 4992 h |
| RTX 4060 Ti 16 GB | $547 | 3006 | 4067 | **550** | 0.105 | 5225 h |
| RTX 5080 | $1530 | 8412 | 11706 | **550** | 0.245 | 6254 h |
| RTX 3090 Ti | $1398 | 6911 | 9939 | **494** | 0.241 | 5791 h |
| RTX 4090 | $2805 | 8925 | 12695 | **318** | 0.417 | 6728 h |
| RTX 5090 | $4600 | 13168 | 19774 | **286** | 0.446 | 10313 h |
| A100 SXM4 40 GB | unknown | 7448 | 11473 | - | 0.402 | - |
| H200 | unknown | 15362 | 25178 | - | 4.739 | - |

## What it says

**Value goes to the old cards.** The top of the table is RTX 2060 Super, RTX 2070, RTX 2070 Super, at 1882, 1874 and 1789 MCells/s per $100. Every one of them is a Turing or Ampere card that the gaming market has finished with, and they run the same HIP kernels as the new ones: the GPU engine needs Pascal or later and nothing else.

**The fastest card for a budget.** Under $300 that is the RTX 2080 Ti (4342 MCells/s on the horn); under $600 the RTX 3080 Ti (5675). At the other end the RTX 4090 at $2805 and the RTX 5090 at $4600 buy their speed at 318 and 286 MCells/s per $100, a third of what the cheap cards give.

**Renting usually wins.** The break even against Vast.ai runs from 1513 to 10313 hours, 3476 in the middle: that is months of continuous simulation before a purchase pays for itself, and the rented machine is somebody else's problem when it breaks. Buying makes sense for a workstation that simulates every day, or where the data cannot leave the building.

**The memory of the AI market sets these prices, not the speed of the cards.** The 24 GB cards (RTX 3090, 3090 Ti, 4090) and the 16 GB RTX 4060 Ti carry a premium for their memory that this workload does not use: openEMS needed under 850 MiB of GPU memory for every benchmark here, so a cheap 8 GB card runs the same job as a 24 GB one. Buy for bandwidth, not for capacity, until a mesh no longer fits.

## The prices

eBay does not serve its sold listings to a script (HTTP 403 without a session), so every price here is a tracker's aggregation of completed eBay sales, not listings counted one by one: getpcparts.com for the older cards, pcpartvalue.com (medians of completed sales) for the GeForce 40 series, cross-checked against jawa.gg, a used marketplace that publishes its own completed sales, and against active asking prices. Each figure blends the models of a card, reference and premium. The used market was rising fast when this was collected, 5 to 8 percent over 30 days for the older cards and 10 to 27 percent over 90 days for the 40 series, so read the table as a snapshot of September 2026. The GeForce 50 series is a special case: only active asking prices were available for it, in the middle of a memory shortage that has every 50 series card 30 to 130 percent above its launch price, so those rows carry low confidence. Where the used asking price stands above new retail, as it does for the RTX 5090, the table takes the retail price, because that is what the card really costs to get.

- **A100 SXM4 40 GB**: no price could be verified (no completed sale could be verified; the band comes from a buying guide that names no source, and the SXM4 40 GB barely trades (the 80 GB card is the volume part). It also needs an SXM carrier board, so it is not a card one simply buys and plugs in)
- **GTX 1080 Ti**: $155 used, $137-177, medium confidence
- **GTX 1660 Ti**: $103 used, $88-103, medium confidence
- **H200**: no price could be verified (it does trade, but only asking prices are published, from $40k to $78k: no median is honest here)
- **RTX 2060**: $138 used, $119-146, medium confidence
- **RTX 2060 Super**: $168 used, $142-197, medium confidence
- **RTX 2070**: $160 used, $135-201, medium confidence
- **RTX 2070 Super**: $174 used, $148-200, low confidence — the asking prices average $228, a gap of 24 percent against the sold figure, far wider than the 4 to 11 percent of the other cards: the clearing price is probably $175 to $200
- **RTX 2080 Ti**: $271 used, $232-300, medium confidence
- **RTX 3050**: $175 used, $151-190, medium confidence — dearer than the faster RTX 2060 and near the RTX 2070: it is still sold new (about $259), which anchors the used price to retail rather than to performance
- **RTX 3060**: $295 used, $212-295, medium confidence — pcprice.watch gives EUR 249 over 1447 completed sales, which agrees
- **RTX 3060 Ti**: $248 used, $228-269, medium confidence
- **RTX 3070**: $260 used, $208-345, high confidence — the best sampled card of the set, and pcprice.watch agrees over 1222 sales
- **RTX 3070 Ti**: $309 used, $270-317, medium confidence
- **RTX 3080 10 GB**: $358 used, $304-391, high confidence — three trackers land within a few percent: $357.59, $358 and EUR 354
- **RTX 3080 Ti**: $487 used, $419-562, medium confidence — the other tracker had two sales only, so this one is the headline
- **RTX 3090**: $1032 used, $1032-1499, low confidence — the trackers disagree badly for the same month ($1032 to $1499) and one 90 day average sits 40 percent under its own 30 day average; the 24 GB cards trade on the demand for memory, not on their speed
- **RTX 3090 Ti**: $1398 used, $1254-1529, medium confidence — at or above its launch price, like the RTX 3090
- **RTX 4060**: $289 used, $275-304, medium confidence — jawa.gg gives $259 over 12 months
- **RTX 4060 Ti 16 GB**: $547 used, $411-650, low confidence — the sources disagree by a third ($411 jawa.gg, $547 pcpartvalue, $650 asking); the 16 GB card goes for about 65 percent more than the 8 GB one ($330), bought for the memory rather than for games
- **RTX 4070**: $515 used, $475-544, medium confidence
- **RTX 4070 Super**: $550 used, $537-617, medium confidence
- **RTX 4070 Ti**: $688 used, $616-752, medium confidence — two trackers of completed sales agree within $10 ($688 and $698); jawa.gg gives $632 over 12 months, lower because it averages the whole year of a rising market
- **RTX 4070 Ti Super**: $801 used, $735-1020, low confidence — three sales only, and the sources span $735 to $1020: read it as $800 to $900
- **RTX 4080**: $905 used, $859-954, high confidence — the one solid figure of the 40 series: jawa.gg independently gives $892 over 34 sales
- **RTX 4080 Super**: $1095 used, $964-1150, medium confidence
- **RTX 4090**: $2805 used, $2500-3010, medium confidence — about 1.75 times its launch price of $1599, which four trackers agree on: the memory of the AI market
- **RTX 5060**: $370 used, asking, $352-388, low confidence
- **RTX 5060 Ti 16 GB**: $640 used, asking, $560-720, low confidence — the pool of listings mixes the 8 GB and 16 GB cards, so the figure is biased low for the 16 GB one
- **RTX 5070**: $700 used, asking, $668-732, low confidence
- **RTX 5070 Ti**: $1050 used, asking, $1002-1098, low confidence
- **RTX 5080**: $1530 used, asking, $1421-1639, low confidence — completed sales across seven markets give about $1490, which agrees
- **RTX 5090**: $5099 new retail, $4620-5579, low confidence — the used asking price ($5099) is ABOVE new retail ($4600), a scalper premium on a card that is hard to buy: the table takes the retail price
