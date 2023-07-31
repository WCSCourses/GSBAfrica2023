Technical and VM queries

**Q. What are the minimum technical requirements?**

To fully benefit from the virtual format, participants will require
minimum computational

specifications and internet access. You can use a PC running Windows or
Linux, or a Mac running OSX. In all you will need admin rights to the
machine in order to install the VM. The minimum recommended requirements
are:

• Webcam and microphone (for interactive video-conferencing, e.g. Zoom)

• 8GB RAM

• i5 processor

• 80 GB free hard disk space

• 8Mbps download speed

**Q: What is the size of the VM? Can I run it from an external drive?**

A: The VM is \~24 GB zipped, but unpacks to \~60GB unzipped. It can run
on an external HD, but this must be a high speed connection (usb3/usb c)

**Q: Do we download virtual box first? then VM?**

You will have to download the virtual box to launch the VM you download
in line with the OS you are running on your computer/ laptop

**Q: I have just followed the provided tutorial and installed the
virtual box to the very end. But when I start up the virtual machine
(ubuntu) I am only seeing a black screen. What could be the issue?**

![](media/image4.png){width="4.433994969378828in"
height="2.507517497812773in"}

A: This needs to be dealt with on a case by case basis as there might be
a few reasons and solutions.

-   Not enough RAM allocated.

-   Possibility of a corrupted download

-   Try increasing the guest\'s video ram

-   Try rebooting in headless mode and once booted select \"Show\"

-   Otherwise send your VMs log that will be somewhere here:
    > \~/VirtualBox VMs/ngs2021/Logs/VBox.log . We can then have a look
    > and see if there is anything we can identify as the problem.

**Q: We have serious problems downloading the VM. For someone who
already uses Ubuntu 64-bits, is it necessary to download the VM?**

A: Yes. There is a lot of software and data on the VM which they won\'t
have.

**Q: Is there a link or a resource on Vula to download the virtual
machine? or do we just download any preferred virtual machine online?**

A: There are announcements on Vula giving details of how to download the
virtual machine. Please do not just use any virtual machine as it won't
have the software and data required for the course.

**Q: My PC is running windows, and I have successfully downloaded the
file through the link, but I can\'t install the extracted file named
NGSBioRC_2021.vdi. Please how do I progress from here?**

**My analysis gave that also: The VM was set up with virtual box version
6.0.x which is no longer supported. So I select a higher version 6.1.
But when I have to choose the system type, I cannot select the Ubuntu
(64-bits) only ubuntu 32-bits was available. Also MD5 number didn\'t
match**

![](media/image2.png){width="5.761111111111111in"
height="3.2402777777777776in"}

A: With regard to the md5, did you run the check against the zipped or
unzipped file?

A: The issue with 32 bit VM is detailed in the document about common
issues you may see during installation. Essentially, you need to go to
your BIOS and change a setting.

**Q: The file is too large and it's taking a lot of time\...**

A: It needs to be that large to include all of the data and software
needed for the course. Be patient, it may take several hours to download

**Q: [Do linux users also required to install the VM?]{.mark}**

[A: Yes. There is a lot of software and data on the VM which they won\'t
have.]{.mark}

**Q: My Virtual Machine is not loading on the Virtuabox platform. I have
tried all thing possible as stated on the documentation provided with
the download. I have another virtuabox machine and it loads. Is it
possible to carry on with that one?**

A: No. The virtual machine for the course contains required software and
data which won't be available on your virtual machine.

**Q: I have been working on Linux using another VM (CentOS 8 64-bit)
that I power it on using program called (VMware workstation pro). is it
possible to use this machine through the course or do I need to download
the VM you mentioned ?**

A: You will need to download the mentioned VM. Your VM will not contain
the course data or all of the software required for the course.

**Q: How can I solve the problem below**

![](media/image1.png){width="5.101253280839895in"
height="3.900161854768154in"} I

A: That looks like you don\'t have admin rights on your PC. Either
install as the admin user or ask someone with admin rights to the
install.

**Q: Can a Staff download the VM onto a pendrive for the participants to
copy instead of each person downloading on their own? This will also
minimize band width challenges I believe**

A: This is fine to do. I wouldn\'t advise running from the pen drive
though.

**Q: I downloaded the vm file but I can\'t open them neither on ubuntu
nor windows. as I couldn\'t download the oracle and I dunno the command
line for unubtu.**

A: You will need to install VirtualBox even on Ubuntu.

**Q: Can i download and install it on external hard drive and work from
it??**

A: You can. The external drive will need to be USB3.0 or better and I
wouldn\'t recommend running from a pen drive.

**Q: How much hard drive space is needed? I have a 64-bit OS, with 4GB
RAM,. How much of hard drive space will be needed for the
installation?**

A: I think 60 to 80GB will be enough space. The VM file initially
expands out to 38GB but will grow as the course goes on.

**Q: Do linux users also required to install the VM?**

A: Yes. There is a lot of software and data on the VM which they won\'t
have.

**Q: I am not able to see working VM. I got message Auto Capture keyboar
optons turned on. This wil cause the Virtual Machine to automatically
capture the keyboar every time the VM winow is activated and make it
unavialable to other applications running on your host machines: When
the keyboar is captured, all keystrokes including system ones like
Alt-TAB) will be directed to the VM. Abd U see image like the
following.**

![](media/image3.png){width="4.62650699912511in"
height="2.9650535870516186in"}

A: The message about autocapture isn\'t something to worry about. You
can close that. However it looks like your VM has booted into grub. This
shouldn\'t happen automatically. Please try closing it and restarting as
an initial step.

**Q: Can I use my linux Ubuntu 2020 machine directly? I have make a
directory of NGS-2021 in my machine. So I don\'t know if it is necessary
installing and working on VM.**

A: Unfortunately, you will need to install VirtualBox and download the
virtual machine file. It contains software and data which you won\'t
have on a standard Ubuntu 2020 LTS build.

**Q: Some colleagues with 64 bit processor laptop were unable to load
Ubuntu64 bit during installation. Please what could be responsible for
that.**

A: I got the same issue, then I download ubuntu it helped

**Q: So if we have already downloaded virtual machine in other projects,
we will have to download it again for this course?**

A: Yes, you have to download the VM specific for this course. It has all
the software and materials for this course.

**Q: HP Workstation Issue**

If anyone is running the VM on an HP Workstation running Windows 10,
they may see the following error:

"Error relaunching VirtualBox VM process: 5"

Followed by:

 "Please try reinstalling VirtualBox.

where: supR3HardenedWinReSpawn what: 5 VERR_INVALID_NAME (-104) -

Invalid (malformed) file/path name. \[/i\]"

It can be resolved by uninstalling the HP Sure Click security software,
followed by restarting the workstation.

**Q: Problem with new Macs**

A: Currently, the main software we use for this course, VirtualBox, is
not compatible with the new Apple M1 chipset. If your computer is a new
Mac using this chipset, please try to source a machine using an Intel
processor.

**Q: VM image works only on 64 bit processor. Unable to have it working
on 32 bit processor on window 7**

The downloaded Ubuntu image will only work on a 64 bit processor. Please
check if virtualisation is enabled in your BIOS and try again.  

 

**Q: I have issues with the file downloads. I use a Mac OS and can\'t
seen to open any of them as expected even with unarchiver. Please is
there a .dmg for mac users OR anything else?**

You need to install and open VirtualBox first. After that go through the
instructions you have downloaded on how to attach the NGSBioRC_2021.vdi
image.

 

**Q: After downloading the VirtualBox through the following
link: [[https://www.virtualbox.org/wiki/Download_Old_Builds_6_0
\[virtualbox.org\]]{.underline}](https://urldefense.proofpoint.com/v2/url?u=https-3A__www.virtualbox.org_wiki_Download-5FOld-5FBuilds-5F6-5F0&d=DwMFaQ&c=D7ByGjS34AllFgecYw0iC6Zq7qlm8uclZFI0SqQnqBo&r=fZ-afGg8pSU0q5yUmvoj9gwq7hVkBMl1Hkj2DAimQOff60ykDPlZHWDIoCETJagh&m=rwiE5ls8QGQuQSN25escldJGh8ZaiVtjBKSSNyeYwKY&s=vHnauEZLMAA0AVd26C1m9XYHv-os9IAaKi5L1dP7eV8&e=) The
next step is downloading the installer for the current operating system.
Which installer is meant?**

On the page there is an option to download an installer for every OS.
Just select the appropriate installer according to your OS.

**VM Download queries**

**Q: During the process of downloading and after 4G download, the
process was cancelled and I get a failed message with no file**

**Q: Every time you download VM and it stops downloading, does it start
from zero when you say resume**

**Q: I tried to download the VM. But at the end downloading, I received
a \"Fail - Forbidden notification.**

If it is related to Google then here are some
suggestions: [[https://appuals.com/how-to-fix-failed-forbidden-error-when-downloading-from-google-drive/
\[appuals.com\]]{.underline}](https://urldefense.proofpoint.com/v2/url?u=https-3A__appuals.com_how-2Dto-2Dfix-2Dfailed-2Dforbidden-2Derror-2Dwhen-2Ddownloading-2Dfrom-2Dgoogle-2Ddrive_&d=DwMFaQ&c=D7ByGjS34AllFgecYw0iC6Zq7qlm8uclZFI0SqQnqBo&r=fZ-afGg8pSU0q5yUmvoj9gwq7hVkBMl1Hkj2DAimQOff60ykDPlZHWDIoCETJagh&m=-Y0T64FTNBhakeoT9I0ioPWTEAtk-xVK8PLXxn8X9HA&s=lsV3K1_WvIR5U_OHlUqL9hQOZyAOZp9Xidpz2Mo2DSQ&e=)

Try the Globus method if direct download links are not working.
