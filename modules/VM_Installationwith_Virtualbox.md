Installing VirtualBox and Setting Up the Virtual Machine for the Viral
Genomics and Bioinformatics 2022 course

**Introduction to VirtualBox and the Virtual Machine**

VirtualBox is freely available software which allows a user to run a
virtual machine on their host computer, meaning that they can run, as we
will, a Linux operating system on Windows or OSX. The virtual machine
contains all of the software and data to be used during the course.
There is also a video of this tutorial, please take note of the name of
your virtual machine (ViralBio22_LAC.vdi.gz).

**Installation of VirtualBox**

VirtualBox is available from this website:

[[https://www.virtualbox.org/]{.underline}](https://www.virtualbox.org/)

Please download the installer for the operating system of your computer
and run through the installation process. **PLEASE NOTE:** admin rights
to your computer are required to install VirtualBox.

Also, **please note** that VirtualBox will not, as yet, run on a new
Apple Mac with the M1 chip set.

After installation of VirtualBox, you should also install the VirtualBox
extension pack, which is available on the same page. The installation
process for this should start after download. The extension pack is not
essential but does include many useful features e.g. the ability to
attach a USB drive to the virtual machine.

Once VirtualBox and its extension pack are installed, running the
software should give a manager window similar to this:

![](media/image8.png){width="6.263888888888889in"
height="4.697916666666667in"}

Note that this is the OSX version of VirtualBox. The Windows interface
differs slightly.

**Download of the Virtual Machine**

The virtual machine is stored as one large file. It can be downloaded
using Globus, and a guide for this is in the folder.

There is also a video (please note, the name of your course VM is:
ViralBio22_LAC.vdi.gz)

The md5 checksum code for verifying that the VM has downloaded correctly
is:

b26528816f3113635445b26e0c8945e8

**Setting Up the Virtual Machine**

PLEASE NOTE: The example images below are taken from previous/different
courses.

Once the virtual machine file has downloaded successfully, you can start
to set the virtual machine up. To do so:

1\) Click the New button as seen on the above screenshot.

2\) This leads to a window similar to this:

![](media/image7.png){width="6.263888888888889in"
height="4.697916666666667in"}

Give a name to your machine. This can be anything but it makes sense to
name it with something related to the course such as 'Viral Bioinf2022.

Next, set the Type to 'Linux' and the version to 'Ubuntu (64-bit)'.
Ubuntu is the precise version of Linux we use for this course. Once
done, click the Continue button.

3\) This leads to a window allowing you to set the memory (RAM)
allocated to the virtual machine:

![](media/image5.png){width="6.263888888888889in"
height="4.697916666666667in"}

The amount of memory you may allocate differs from machine to machine,
depending on how much your computer has. For the course, it's best to
set the memory as high as possible whilst leaving enough to run your
computer's operating system. Thus, we advise that you set the memory so
that it is close to the top end of the green part of the line in the
above screenshot but not so high that the pointer reaches the pink part
of the line. Once the memory is set, click the Continue button.

4\) The next thing to do is to point VirtualBox to the virtual machine
file you have downloaded. This is done using this Window:

![](media/image1.png){width="6.263888888888889in"
height="4.697916666666667in"}

Select the 'Use an existing virtual hard disk file' and click the yellow
and green icon to the right of the pulldown menu. This should lead to a
window similar to this:

Click the Add button on the top left and choose the .vdi file you have
downloaded. This file will then appear in a list of Unattached files in
the above window. Highlight it and click Choose.

The next step after this is to click the create button. If all has gone
correctly, the manager window should now show looking like this:

![](media/image2.png){width="6.263888888888889in"
height="4.697916666666667in"}

5\) As a last step in setting up the virtual machine for the course, we
need to set the number of processors it can access. To do this,
highlight your VM in the manager window as above and click the Setting
button. Choose System then Processor. Similar to setting memory, set the
number of processors (CPUs) used to the upper limit of the green part of
the line:

![](media/image3.png){width="6.263888888888889in"
height="4.942361111111111in"}

As with memory, the number of processors you may use will depend on your
host machine.

You may also find that changing the video memory under the display tab
to 128MB will help avoid black screen and other display issues.

![](media/image4.png){width="3.468144138232721in"
height="3.2656255468066493in"}

6\) All that remains now is to start the virtual machine. On the manager
window, highlight your virtual machine's name and click the 'Start'
icon. The virtual machine should run through a boot process and, after a
short time, you should see a window similar to this:

![](media/image6.png){width="6.263888888888889in"
height="3.640277777777778in"}

The user account on the virtual machine is named 'manager' and the
password, if required, is also 'manager'.

Check that the key board has mapped correctly to the virtual machine

/ \\ \| @ \~ & \#

In the VM, check the settings, then keyboard settings, and check that
the keyboard versions and locations are matched.
