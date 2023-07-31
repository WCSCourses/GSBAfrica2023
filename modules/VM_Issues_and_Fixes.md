---
title: "[]{#_gynl14krzeyx .anchor}Summing up issues and fixes on viral
  bioinf 2022 VM:"
---

Data on phylogeny and typing has been changed

[[https://github.com/WCSCourses/ViralBioinfAsia2022/tree/main/course_data/phylogeny_typing]{.underline}](https://github.com/WCSCourses/ViralBioinfAsia2022/tree/main/course_data/phylogeny_typing) -
In the earlier version in one of the input fasta files the \'time of
isolation\' information was missing.

**To do:** pull down new data onto VM with git pull and lfs pull for
final VM

From richard on sarscov-2 workflow 1-08-2022

Had a number of issues when I was testing things and got some
workarounds (nothing wrong with the VM, more finding bugs in all the
tools, some of them quite major ones!)

-   Minor point. The data is currently duplicated. There is a full
    > unzipped version of the data in \~/SARS-CoV-2/, and then a still
    > zipped one in \~/course_data/SARS-CoV-2_workflows/SARS-CoV-2. I
    > used the data in \~/SARS-CoV-2 for testing as already unzipped. So
    > if you want to save space you could delete
    > \~/course_data/SARS-CoV-2_workflows/SARS-CoV-2 and/or move
    > \~/SARS-CoV-2 to \~/course_data/SARS-CoV-2_workflows/SARS-CoV-2 --
    > let me know which one is going to be used so I can get the turoial
    > in the right folder.

Matt Basthon,

FYI I don't have \~/course_data/ in my VM at all, not sure if the image
has changed since Richard downloaded it?

-   SPEAR -- didn't work, crashed with an error all the time, ended up
    > contacting the developer and turned out to be a major but simple
    > bug and they fixed it and updated github, I managed to update my
    > VM by simply running from within the spear conda env (so conda
    > activate spear first) '\~/SPEAR/scripts/update_spear.sh spear'
    > github said to do spear update spear but that didn't work for me:

```{=html}
<!-- -->
```
-   <https://github.com/m-crown/SPEAR>

Matt Bashton, run

/home/manager/miniconda/envs/spear/bin/update_spear.sh spear

Which will update SPEAR, it's not clear why spear can't find the
update_spear.sh script as it's in the path.

-   Civet -- again didn't work when testing (it's all installed fine) --
    > it looks to me that the the tool gofasta from the same Edinburgh
    > group is now incompatible with Civet (and has been for some time I
    > think) I won't bore with details. I put an issue on GitHub, but
    > the workaround I did to get things to work correctly was to
    > download an older version of gofasta and replace the current
    > version in the civet conda env.

Matt Bashton: If you pin the older version into the conda env then it
should install only the older version, not sure of version numbers, but
there was a lot of pre 1.0 verisons.

-   I downloaded the "13th July 2021" version from gofasta github
    > gofasta-linux-amd64:
    > <https://github.com/virus-evolution/gofasta/releases>

    > And then copied that file to replace the one in the civet conda
    > env (make a backup first though):

    > cp \~/miniconda/envs/civet/bin/gofasta \~/

    > cp \~/Downloads/gofasta-linux-amd64
    > \~/miniconda/envs/civet/bin/gofasta

```{=html}
<!-- -->
```
-   NextFlow viralrecon -- a complete mare from start to finish. Had
    > numerous issues, sorted many of them out (mostly how to run it,
    > it's defaults are looking for 6 cpus and 36GB of RAM and all this
    > over stuff). Even after sorting that it doesn't work properly,
    > hangs alot. So going to ditch it, it's an overkill and probably
    > too ambitious for a 1.5 hour session to introduce properly. So
    > going to go back to the original illumina plan (pre-netflow) of
    > teaching the four main step trim_galore/bwa/ivar trim/ivar
    > consensus -- I tested them and they work fine -- this is probably
    > more interesting for people to learn anyway rather than an opaque
    > NetFlow pipeline.

```{=html}
<!-- -->
```
-   But I am going to see if I can get the COG-UK sars nextflow working
    > easily (the one with the DSL error) as this is simple and easy to
    > use and just uses the same steps as I was going to teach anyway. I
    > not basing the course around it, but if I can fihure out how to
    > get it to work (I'll contact the developer) would be good to have
    > it working on the vm (i.e. a might be able to write a script to
    > students to fix it later on).

    > I'll point people to viralrecon and the cambridgebio tutorial as a
    > way for doing 100's of samples (which is prob out of scope of the
    > VM anyway)

The artic-ncov19 conda env worked completely fine. (edited)

Matt Bashton, Let me know if you have any issues with the Connar lab
nextflow, I've used it a lot over the last two years, so know it's
inside out if it's not playing ball. Sometimes there are issues with
versions of stuff in conda breaking the install but you can normally get
round that by pinning the versions.

[[Martin Aslett]{.underline}](mailto:maa@sanger.ac.uk) October has added
data into Reference alignment and Consensus n Variant calling course
data, pls lfs pull the git again on VM

I see David is using -
[[https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/Modules/Coverage_plots_and_statistics.md]{.underline}](https://github.com/WCSCourses/ViralBioinfAsia2022/blob/main/Modules/Coverage_plots_and_statistics.md)
this breadcrumb - cd
/home/manager/ViralBioinfAsia2022/course_data/Coverage_Plots_Stats

So lets keep course data as is.
