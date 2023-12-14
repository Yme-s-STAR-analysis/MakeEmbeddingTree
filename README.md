# Efficiency Step 1

Get embedding reduced tree

`version 3.0`

`author: yghuang`

***

## Quick Start

1. Copy this dir. and rename it (like `proton`, `antiproton`, ...). Change cuts in `StMuAnalysisMaker`.

2. `cons`.

3. Run `./prepare.sh $path`, here `$path` is the target dir. to save the reduced tree.

4. Change the file list in `Csubmit.xml` and `star-submit Csubmit.xml`.

5. Before submitting the jobs in step 4, you can prepare a small demo in `test.list` and directly run `./runtest.sh` for testing.


## Change log

14.12.2023 by yghuang (v3.0):

> StRefMultCorr package will be used in this new patch.
>
>> The Centrality bin will be saved instead of multiplicity.

14.08.2023 by yghuang (v2.1):

> Using new centrality utils.
>
> For different data set, replace `Conf.h` in `CentralityUtils` and `cons`.

22.08.2022 by yghuang (v2.0):

> Reimplemented.

old version

2022 Mar. 30 by yghuang (1.5):

> Reimplemented.

2022 February 10 by yghuang (1.4):

> Update: changed for 19.6GeV use.
>
> Outlook: In the future, may use pico dst instead of mu dst.

2021 September 17 by yghuang (1.3):

> Update: cancel the trigger selection.

2021 September 13 by yghuang (1.2):

> Update: save more event information for cutting, instead of doing the cuts when read MuDst.
>
> > (Which means the codes used in next step should add those cuts.)

2021 September 12 by yghuang (1.1):

> Update: change the method of submitting jobs. `./Ssubmit.sh` -> `./prepare.sh` + `star-submit Csubmit.xml`.
>
> Update: modify the codes in StMuAnalysisMaker, which maintains less useless variables for saving.

2021 September 9 by yghuang (1.0):

> A test version, not finished.