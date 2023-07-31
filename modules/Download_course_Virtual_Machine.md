**Delivery of the course Virtual Machine (VM) image**

The VM is a large download that many users struggle with.

We are going to use Globus software to help with this. Globus is file
transfer software designed for research. Here is a
[[video]{.underline}](https://drive.google.com/file/d/1sr8bRoiq0XS6EJ09UWbu-Ag6MW1LB6IR/view?usp=sharing)
of this tutorial.

Globus website:

[[https://www.globus.org]{.underline}](https://www.globus.org)

VirtualBox is used to run the VM once it is downloaded.

[[https://www.virtualbox.org]{.underline}](https://www.virtualbox.org)

Steps:

1.  Make a globus connect personal account by going this webpage, and
    > choosing the download for your operating system, which will prompt
    > you to make an account:
    > [[https://www.globus.org/globus-connect-personal]{.underline}](https://www.globus.org/globus-connect-personal)

![](media/image6.png){width="3.744792213473316in"
height="2.519767060367454in"}

Choose use "Globus ID to sign in" (on bottom)

![](media/image8.png){width="3.838542213473316in"
height="1.7853685476815397in"}

Then if you dont have an ID yet, select "Need a Globus ID? Sign up"

![](media/image4.png){width="3.0012817147856516in"
height="2.7968755468066493in"}Then make sure you specify for research or
educational purposes. And create your account. Remember your password
for the later steps.

2.  Then download the globus client onto your local machine (or where
    > you intend to run the VM), allow it to install. It will ask for a
    > collection name, give it a name you will refer to - example
    > "home_computer" or "local_mac". This is a name for the local
    > folders on your computer that we will send the VM to in a later
    > step.

3.  Click on the small "g" icon on the task bar and select Web: Transfer
    > Files. For linux users - there may not be a shortcut. Once you
    > start globus personal connect via command line, navigate to
    > [[[https://app.globus.org/file-manager]{.underline}]{.mark}](https://app.globus.org/file-manager)
    > to begin the file manager.\
    > ![](media/image5.png){width="3.0781255468066493in"
    > height="0.9922900262467191in"}

> ![](media/image3.png){width="2.0088101487314085in"
> height="2.2918700787401574in"}

4.  You will see your own files on this page. Your local endpoint is
    > your computer. Click on "Collections" on the left

> ![](media/image2.png){width="1.5160870516185476in"
> height="3.0677088801399823in"}

5.  Search for the endpoint at **wcs_data_transfer**

> ![](media/image10.png){width="6.267716535433071in"
> height="1.3888888888888888in"}
>
> Click on the endpoint labelled **wcs_data_transfer**
>
> ![](media/image9.png){width="4.744792213473316in"
> height="1.5553466754155731in"}
>
> Then click on "Open in File Manager" will open up a file manager

6.  Then we begin the steps to transfer the VM to your local machine.
    > First, select the VM file "ViralBio22_LAC.vdi.gz" with the check
    > box. Then click on "Transfer or sync to"

> ![](media/image1.png){width="6.267716535433071in"
> height="2.5277777777777777in"}
>
> Then click on the search box in the opposite panel
>
> ![](media/image7.png){width="6.267716535433071in"
> height="1.5694444444444444in"}
>
> Click on the "Your Collections", and select your local endpoint (the
> name will be what you gave it during the globus personal connect
> installation)
>
> You can browse to choose the specific directory or folder on your
> local machine you want the VM to download to.
>
> ![](media/image11.png){width="6.267716535433071in"
> height="3.013888888888889in"}

7.  When you have chosen the local location, Click on the "Start" button
    > under the "**wcs_data_transfer**" section to begin the download to
    > your local endpoint.

8.  Wait for download completion - it will email you to the account you
    > set up, and you can track the transfer in the "Activity" menu.

9.  Run the installation of virtualbox, then install the virtual machine
    > you have just downloaded. See install docs on the drive, including
    > VM install video.

Note - if you encounter severe problems with globus, such as it being
blocked on an institutional firewall, then there is a backup option to
download via google drive:

**\
\
**

**There is a google drive link to download the VM directly:**

**This for the Vbio LAC 2022 VM:**

**[[https://drive.google.com/file/d/1GB7MKikyKp82OFqC-IJdC_vuRoUryn9N/view?usp=sharing]{.underline}](https://drive.google.com/file/d/1GB7MKikyKp82OFqC-IJdC_vuRoUryn9N/view?usp=sharing)**

NB Caution, direct downloads can break easily, and lead to corrupted
VMs, use the md5 checksum to check if your VM is as it should be.

to do an md5 check -
<https://portal.nutanix.com/page/documents/kbs/details?targetId=kA07V000000LWYqSAO>

the md5 must be: b26528816f3113635445b26e0c8945e8
