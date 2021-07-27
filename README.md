PanGIA GUI Manual
================
7/14/2021

<figure>
<img src="gui/images_gui/pangia_logo.png">
<figcaption>
<em>Your Guide to the PanGIA GUI</em>
</figcaption>
</figure>

## **Table of Contents**

### 1. [Introduction](#section-1-introduction)

### 2. [Initial Download](#section-2-initial-download)

      **2.1**: [**PanGIA & PanGIA
GUI**](#section-21-pangia-and-pangia-gui)  
          *2.1.1*: [*Downloading from the Command
Line*](#section-211-downloading-from-the-command-line)  
          *2.1.2*: [*Downloading
Manually*](#section-212-downloading-manually)  
      **2.2**: [**Reference
Databases**](#section-22-reference-databases)

### 3. [Preparing the Compute Environment](#section-3-preparing-the-compute-environment)

      **3.1**: [**Preparing a Non-Dockerized
Environment**](#section-31-preparing-a-non-dockerized-environment)  
      **3.2**: [**Preparing a Dockerized
Environment**](#section-32-preparing-a-dockerized-environment)

### 4. [Using the GUI](#section-4-using-the-gui)

      **4.1**: [**Accessing the GUI**](#section-41-accessing-the-gui)  
          *4.1.1*: [*Non-Dockerized
Case*](#section-411-non-dockerized-case)  
          *4.1.2*: [*Dockerized*](#section-412-dockerized-case)  
      **4.2**: [**GUI Navigation**](#section-42-gui-navigation)  
          *4.2.1*: [*Dashboard*](#section-421-dashboard)  
          *4.2.2*: [*Users*](#section-422-users)  
          *4.2.3*: [*Settings*](#section-423-settings)  
              *4.2.3.1*: [*General
Settings*](#section-4231-general-settings)  
              *4.2.3.2*: [*Project
Settings*](#section-4232-project-settings)  
              *4.2.3.3*: [*File Template
Settings*](#section-4233-file-template-settings)  
              *4.2.3.4*: [*Meta Type
Settings*](#section-4234-meta-type-settings)  
          *4.2.4*: [*Projects*](#section-424-projects)  
              *4.2.4.1*: [*Running PanGIA from the Projects
Tab*](#section-4241-running-pangia-from-the-projects-tab)  
              *4.2.4.2*: [*Exploring the Projects
Tab*](#section-4242-exploring-the-projects-tab)  
          *4.2.5*: [*PanGIA and PanGIA
Results*](#section-425-pangia-and-pangia-results)  
              *4.2.5.1*: [*Results Page*](#section-4251-results-page)  
          *4.2.6*: [*Real-Time PanGIA*](#section-426-real-time-pangia)  
      **4.3**: [**Examples of Usage**](#section-43-examples-of-usage)  
          *4.3.1*: [*Baseline Example*](#section-431-baseline-example)  
          *4.3.2*: \[*Example \#2: Real-Time Case with T-MARK
analysis*\]  
          *4.3.3*: \[*Example \#3: Decision Tree - All Results Mode*\]

### 5. [Planned Visualizer Features](#section-5-planned-visualizer-features)

### 6. [Troubleshooting](#section-6-troubleshooting)

      **6.1**: [**Installation
Troubleshooting**](#section-61-installation-troubleshooting)

## *Section 1)*: **Introduction**

<div style="text-align: justify">

This report contains instructions for the installation and use of the
graphical interface designed for PanGIA. PanGIA is standalone and may be
run independently from the command line. However, this forces the
responsibilities of manual configuration and record keeping for each run
entirely onto the user.

This GUI companion to PanGIA makes the pathway to discovery more
efficient and exciting by streamlining the run specification/execution
process, providing a robust job queue, enabling real-time PanGIA
processing, and constructing an intuitive, informative, dynamic
visualization tool for result analysis. The GUI is composed of six
parent pages: Dashboard, Users, Settings, Projects, PanGIA, and Real
Time. If you have already completed PanGIA setup and are using this
document as a GUI reference, please consult [**Section
4.2**](#section-42-gui-navigation) onward.

</div>

## *Section 2)*: **Initial Download**

This section details the process of downloading PanGIA, the GUI, and
curated reference databases necessary to run PanGIA.

### **Section 2.1)**: **PanGIA and PanGIA GUI**

#### *Section 2.1.1)*: *Downloading from the Command Line*

<div style="text-align: justify">

1.  In a terminal, navigate to preferred local path.

2.  Git clone the repo from
    <https://github.com/mriglobal/pangia_gui.git> - run the command:
    `conda install -c anaconda git` if the git package is not present.

#### *Section 2.1.2)*: *Downloading Manually*

1.  Go to the Github page: <https://github.com/mriglobal/pangia_gui>

2.  Download the repository and extract the zipped files to your
    preferred local path.

### **Section 2.2)**: **Reference Databases**

**Talk about ftp here**

## *Section 3)*: **Preparing the Compute Environment**

PanGIA requires a specific compute environment to satisfy its run-time
dependencies. We consider two cases: an environment constructed locally
within a specific Linux distribution, and a PanGIA-friendly Docker
environment independent of the host machine.

### **Section 3.1)**: **Preparing a Non-Dockerized Environment**

1.  Install conda and the correct Linux distribution. See XXXX for
    details.
    1.  We recommend that you use the *Miniconda* distribution available
        at <https://docs.conda.io/en/latest/miniconda.html>. If you’d
        like to download the distribution directly, please use this
        link:
        <https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh>
2.  In a terminal, navigate to PanGIA GUI directory created in
    [**Section 2.1**](#section-21-pangia-and-pangia-gui). Please note
    that the default relative path is ***`pangia_gui/gui`***. Use the
    command: `conda env create –f  gui.yml`. This command acts on the
    .yml file within the directory to construct a python environment
    that satisfies PanGIA’s dependencies. By default, name of the
    environment will be ‘gui’, but this can be replaced by renaming the
    prefix of the .yml file prior to creating the environment.

~~3. Within the terminal from the previous step, use the command:
`conda activate XXX`, where XXX is the name of the environment (‘gui’ by
default). Next, run the following series of commands to finish initial
database and environment setup: + `export FLASK_APP=pangia_gui.py` +
`flask db init` + `flask db migrate –m` + `flask db upgrade` +
`conda install -c bioconda fastp`~~

Whenever accessing the GUI, begin by opening ***three separate command
line terminals***, and navigate each into the PanGIA GUI directory.
Next, execute the `conda activate XXX` command, where ‘XXX’ is the name
of the conda environment specified above in [**Section
3.1**](#section-31-preparing-a-non-dockerized-environment) (‘gui’ by
default) in each of the three separate command line instances/windows.
Each instance is responsible for running one of three services
comprising the GUI. See [**Section 4.1**](#section-41-accessing-the-gui)
for more information. The GUI can be completely reset by navigating to
<http://localhost:5000/auth/reset_site> - this is not recommended.

### **Section 3.2)**: **Preparing a Dockerized Environment**

1.  Download Docker and Docker Compose. Compose may need to be installed
    separately, depending on user OS. More details available at XXXX.

2.  On the command line, navigate to the PanGIA GUI directory.

3.  Open `docker-compose.yml` with a command line editor. Alternatively,
    open the file with a text editor.

4.  Notice the spacing: find the first line under services –&gt; app
    –&gt; volumes. Change this line to specify the path of your PanGIA
    directory. This will grant the Docker container access to everything
    within that directory and mount it into the ‘app’ container.

5.  On the command line and within the PanGIA GUI directory, use
    command: `docker-compose build`.

6.  Verify image construction with Docker Desktop or through Docker
    command line tools. More information on preparing images and
    containers can be found at XXXX.

## *Section 4)*: **Using the GUI**

Running PanGIA through the GUI assumes that three independent services –
a worker, scheduler, and Flask server – are already operational.
[**Section 4.1**](#section-41-accessing-the-gui) covers the process of
bringing each service online and opening the GUI. It also considers the
Dockerized case, in which much of this has already been handled in
[**Section 3.2**](#section-32-preparing-a-dockerized-environment) by
Docker Compose. A complete tour of the GUI can be found in [**Section
4.2**](#section-42-gui-navigation). A series of examples for three
typical PanGIA runs can be found in [**Section
4.3**](#section-43-examples-of-usage), including standard analysis
(Section 4.3.1); T-MARK & Decision Tree (Section 4.3.2); and Real-Time;
(Section 4.3.3).

### **Section 4.1)**: **Accessing the GUI**

#### *Section 4.1.1)*: *Non-Dockerized Case*

1.  Following instructions listed at the the end of [**Section
    3.1**](#section-31-preparing-a-non-dockerized-environment), prepare
    three command-line terminals as described.

</div>

2.  In first window use command: `redis-server`.

<figure>
<img src="gui/images_gui/GUI_Setup/redis.png">
<figcaption>
<em>the redis-server in a terminal window</em>
</figcaption>
</figure>

<br>

3.  In second window use command: `rq worker pangia-tasks`.

<figure>
<img src="gui/images_gui/GUI_Setup/Terminal_RQ-Worker.png">
<figcaption>
<em>the worker in a terminal window</em>
</figcaption>
</figure>

<br>

4.  In third window use commands: `export FLASK_APP=pangia_gui.py` +
    `flask run`.

<figure>
<img src="gui/images_gui/GUI_Setup/Terminal_Flask_App.png">
<figcaption>
<em>the flask app in a terminal window</em>
</figcaption>
</figure>

<br>

5.  In a web browser, navigate to localhost:5000 for GUI access.

<div style="text-align: justify">

#### *Section 4.1.2)*: *Dockerized Case*

1.  In any command line window, navigate to PanGIA GUI directory.

2.  Use command: `docker-compose up`.

3.  In a web browser, navigate to localhost:8000 for GUI access.

### **Section 4.2)**: **GUI Navigation**

#### *Section 4.2.1)*: *Dashboard*

After logging in, users are directed to this page. It lists the history
of both the most recent PanGIA runs and ongoing runs. Links to results
and/or the ongoing job status page is provided under each listing. At
the very bottom of the page, users may click a link directing to the
PanGIA section of the GUI, described in more detail below.

</div>

<figure>
<img src="gui/images_gui/Users_and_Dashboard/PanGIA_Dashboard.png">
<figcaption>
<em>the dashboard</em>
</figcaption>
</figure>

<div style="text-align: justify">

### **Section 4.2.2)**: **Users**

Interface for adding users to the GUI’s persisted database. Current
users are listed with name, username, email, date of last login, and
role (Admin, User, etc.).

</div>

<figure>
<img src="gui/images_gui/Users_and_Dashboard/PanGIA_Users.png">
<figcaption>
<em>the users page</em>
</figcaption>
</figure>

<div style="text-align: justify">

### **Section 4.2.3)**: **Settings**

Multi-tabbed interface for specification of general settings, PanGIA
preprocessing/run parameters, project categories, file templates, and
meta-types. Specifics for each tab and sub-tab are listed below.

##### *Section 4.2.3.1)*: General Settings

These settings are dedicated to tailoring PanGIA runs to the user’s
default needs. Most general settings may be adjusted on-the-fly, while
preparing any PanGIA run. However, preparing defaults here will save
time. The general settings are subdivided into three categories:

1.  Application Settings: User specifies PanGIA directory path and
    PanGIA database/uploads path. On Docker, the pathway must reference
    the PanGIA directory mounted to the ‘app’ container. Important: this
    relative path should match the updated line in the
    docker-compose.yml file, as described in [**Section
    3.2**](#section-32-preparing-a-dockerized-environment) above. Users
    may adjust thread-count made available for PanGIA. An even number is
    recommended. Default is two threads.

</div>

<figure>
<img src="gui/images_gui/Settings/PanGIA_Settings_General_AppSettings.png">
<figcaption>
<em>application settings - within general settings</em>
</figcaption>
</figure>

<div style="text-align: justify">

2.  Preprocessing: Options to alter/filter reads prior to PanGIA
    analysis. May trim reads for quality level, adjust the average
    quality, minimum read length, and “N” base cutoff values, apply a
    low complexity filter, exclude PolyA tails, or indicate a custom
    number of base pairs to cut from 5’ or 3’ ends of reads. The user
    may opt to skip preprocessing.

</div>

<figure>
<img src="gui/images_gui/Settings/PanGIA_Settings_General_Proprocessing.png">
<figcaption>
<em>preprocessing settings - within general settings</em>
</figcaption>
</figure>

<div style="text-align: justify">

3.  Run Parameters: Options determining read scoring, and additional
    forms of analysis (TMARK, Decision Tree, etc.) May specify values
    for seed length, and assign minimal values for score, aligned score,
    read count, read RSNB, linear length, percent genome coverage, and
    depth/RS depth. Scoring method may be selected from a drop-down menu
    (Standalone/Background/Combined), and a series of check boxes
    determine if additional analyses are included in the run.

</div>

<figure>
<img src="gui/images_gui/Settings/PanGIA_Settings_Run_Param.png">
<figcaption>
<em>run parameter settings - within general settings</em>
</figcaption>
</figure>

<div style="text-align: justify">

##### *Section 4.2.3.2)*: *Project Settings*

Used to consolidate and organize related PanGIA uploads/runs. Existing
Projects may be archived or deleted. Projects are instantiated with a
name, description, and slug. They may be constructed as children of an
existing parent Project. When providing a name and slug, use only
alphanumeric and ‘\_’ characters.

</div>

<figure>
<img src="gui/images_gui/Settings/PanGIA_Settings_Projects.png">
<figcaption>
<em>project settings</em>
</figcaption>
</figure>

<div style="text-align: justify">

##### *Section 4.2.3.3)*: *File Template Settings*

User may specify PanGIA input file types other than .fastq on this tab.
File Templates require a name, description, and file extension. When
creating a new File Template, the user must decide which Meta-Types will
be available for runs using that File Template. A default for .fastq
files is included.

</div>

<figure>
<img src="gui/images_gui/Settings/PanGIA_Settings_File_Templates.png">
<figcaption>
<em>file template settings</em>
</figcaption>
</figure>

<div style="text-align: justify">

##### *Section 4.2.3.4)*: *Meta-Type Settings*

User may impose additional rules regarding meta-data, grouping it
together and limiting input variance. Default Meta-Types include names
of spiked-in organisms and spike-concentration, sequencing protocol, and
the sampling date. User may specify whether data was a control sample.
New Meta-Types require a name, description, and value (string, integer,
Boolean, date, etc.), and may be designated as ‘required’ if desired.

</div>

<figure>
<img src="gui/images_gui/Settings/PanGIA_Settings_MetaTypes.png">
<figcaption>
<em>meta-type settings</em>
</figcaption>
</figure>

<div style="text-align: justify">

#### **Section 4.2.4)**: **Projects**

Provides the name of each Project constructed in the Settings - Projects
tab as described above in [**Section
4.2.3.2**](#section-4232-project-settings). Clicking on the name of a
Project yields a list of all the files that have been uploaded to that
Project. The list is grouped by File Template. Files may be ordered by
ascending/descending or date/name. Each File Template header has a large
blue button labeled ‘Add New File’ - this directs the user to the upload
interface. Files uploaded through the GUI are copied into a new
directory named after the Project in the local pangia/upload directory.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/PanGIA_Projects_Main.png">
<figcaption>
<em>file list within a project</em>
</figcaption>
</figure>

<br>

<figure>
<img src="gui/images_gui/Projects_and_Runs/upload_files_browser.png">
<figcaption>
<em>file browser - default upload path corresponds to general settings
specification</em>
</figcaption>
</figure>

<div style="text-align: justify">

##### *Section 4.2.4.1)*: *Running PanGIA from the Projects Tab*

To start a job from the Project menu, use the ‘Run PanGIA’ button on any
file uploaded to the Project. Each run requires a name and description.
Drop-down menus specify which paired files are to be used. The first
menu option defaults to whatever file the Run Information page was
accessed from. All General settings may be adjusted here: the interface
is identical to that appearing in [**Section
4.2.3.1**](#section-4231-general-settings).

##### *Section 4.2.4.2)*: Exploring the Projects Tab

Files listed within a Project have three interactive links. The View
Meta Data and View Results links direct the user to a new page, with two
tabs. These tabs are named X Information, where ‘X’ is the File
Template, and Results.

1.  The X Information tab specifies the filename, upload date, filepath,
    and description.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/PanGIA_Projects_XInfo.png">
<figcaption>
<em>the “X” info information tab</em>
</figcaption>
</figure>

<div style="text-align: justify">

2.  The Results tab lists all PanGIA runs that the file has been part
    of. Details include run name/date, description, parent Project,
    user, and filenames.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/PanGIA_Projects_Results.png">
<figcaption>
<em>the results tab within a project</em>
</figcaption>
</figure>

<div style="text-align: justify">

3.  Runs have a link labeled ‘View Results’. Following this link directs
    the user to the PanGIA Results section of the GUI, described below
    in [**Section 4.2.5**](#section-425-pangia-pangia-results).

#### **Section 4.2.5)**: **PanGIA and PanGIA Results**

Lists all queued and completed jobs, regardless of Project. Each
completed job contains the same information found in the Results tab
accessible through the Projects section of the GUI; each job also has
the same ‘View Results’ link. Both are described above in [**Section
4.2.4**](#section-424-projects).

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/PanGIA_Main.png">
<figcaption>
<em>the PanGIA runs page accessed via the GUI sidebar</em>
</figcaption>
</figure>

<div style="text-align: justify">

##### *Section 4.2.5.1)*: *Results Page*

The Results page provides run name and description, and features an
interactive window to browse rows of the .tsv file output of the PanGIA
run. This page also links to the Visualizer, which will graphically
display that output. An in-depth PanGIA run log is provided at the
bottom – the printout there is similar to what can be found in in the
worker printout during the run. An example Results page appears below.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/example_results_page.png">
<figcaption>
<em>the results page of a PanGIA run</em>
</figcaption>
</figure>

<div style="text-align: justify">

#### **Section 4.2.6)**: **PanGIA Real-Time**

User specifies a job name, description, and associated Project, as
usual. Instead of paired files, the user must provide the path to a
directory that Real-Time PanGIA is listening for. This directory should
be the dumping point for ongoing nanopore sequencing. Notably, once
Real-Time PanGIA is initiated, the user must manually terminate the job
as desired. Completed Real-Time jobs persist PanGIA results from the
last job iteration.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/PanGIA_RealTime.png">
<figcaption>
<em>setting up a RealTime PanGIA run</em>
</figcaption>
</figure>

<div style="text-align: justify">

### **Section 4.3)**: **Examples of Usage**

#### **Section 4.3.1)**: **Baseline Example**

The following walkthrough describes a typical PanGIA run with default
settings. For this example, assume the Project `SOP_Runs` has already
been established, and is accessible through the GUI sidebar as shown
below.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/current_project_demo.png">
<figcaption>
<em>project sidebar dropdown menu - leads to project file list</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

This link redirects to the file page associated with ‘SOP\_Runs’. For
the purposes of this tutorial, assume that two .fastq files – named
`BMal_100000_rep1_R1` and `BMal_100000_rep1_R2` – have already been
uploaded to the Project.

Clicking on the Run PanGIA buttons associated with either .fastq file
directs to the Run PanGIA page. As discussed in [**Section
4.2.4.1**](#section-4241-general-settings), this page is split into Run
Information and General settings. Pressing the button associated with
`BMal_100000_rep1_R1` directs to the Run Information page shown below.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/run_info_example_bmal100k_final.png">
<figcaption>
<em>selecting paired fastq files</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

This run has been given a name and description, and a selection for the
`Paired Fastq` has been made – in this case, `BMal_100000_rep1_R2`. None
of the General settings on the Run PanGIA page need to be adjusted, as
this is a default PanGIA run – so this job is ready to be executed by
clicking on the Run PanGIA button at the bottom of the page!

The PanGIA job is now running, which may be verified in a variety of
ways. First, navigate to the command line terminal running the rq
worker. It should look something like this:

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/rq_worker_printout.png">
<figcaption>
<em>typical printout in the rq worker terminal during a PanGIA run</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

The rq worker provides time-stamped logs for run status and allows the
user to follow along with each compute stage, or troubleshoot. Heading
back to the GUI, navigate to the Dashboard. It should look something
like this:

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/ongoing_job_dashboard.png">
<figcaption>
<em>the dashboard dynamically changes while PanGIA runs are
ongoing/queued</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

The current job `BMAL100k` is listed – and the link below the job
directs to the Running Job page. If additional PanGIA runs were queued
behind `BMAL100k`, they would be listed here. The Running Job page is a
run log, and should look similar to the rq worker printout, as seen
below:

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/ongoing_job_job_page.png">
<figcaption>
<em>the job log page - includes similar information to the rq worker
during jobs/queued</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

After the job is finished, the Running Job page turns into a Results
page. Job completion will be reflected in the PanGIA page of the GUI:

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/finished_run_jobs_page_example1.png">
<figcaption>
<em>the jobs page updating following run completion</em>
</figcaption>
</figure>

<br>

<figure>
<img src="gui/images_gui/Projects_and_Runs/updated_finjobs_page_example1.png">
<figcaption>
<em>the jobs page updating following run completion</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

Clicking the PanGIA Visualization button opens the Visualizer in a
separate window. More detail concerning planned features are described
in [**Section 5**](#section-5-planned-visualizer-features), but the
scatterplot for the static PanGIA run as described in this section
should appear similar to the one depicted here:

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/vis_example1.png">
<figcaption>
<em>the visualized results of the run</em>
</figcaption>
</figure>

<div style="text-align: justify">

##### *Section 4.3.1.1)*: *Baseline Example Visualizer Usage and Analysis*

This section discusses current features of the Visualizer and uses the
PanGIA run from [**Section 4.3.1**](#section-431-baseline-example) to
describe how the Visualizer is used to interpret PanGIA results.
Presently, the Visualizer displays hits for taxa between the genus and
strain level, and includes a scatterplot and datatable. On the
scatterplot, each dot on the plot corresponds to one such taxa and are
listed on the x-axis.

1.  The scatterplot y-axis indicates `Raw Read Count` of each hit. This
    axis is dynamically scaled based on maximum/minimum values present
    in the sample. In this example, the y-axis covers a very large range
    of values – between an effective measurement of zero reads and
    one-million. While this would be too large a scale for comparison of
    read count between two strain-level hits, it is a very suitable
    scale in this example, where genus *Burkholderia* is compared
    against species *Burkholderia mallei* and a further strain subset
    *Burkholderia pseudomallei*. In the image above, the center of the
    point representing genus *Burkholderia* intersects the y-axis (read
    count) somewhere between fifty-thousand and five-hundred thousand
    reads – followed in order of magnitude by species and strain.

2.  The scatterplot color-bar `Score-Gradient` axis on the right
    indicates hit performance on the selected `Score` metric – either
    `Standalone`, `Background`, or `Combined`. The color-bar is also
    dynamically scaled; in the above example, strain *Burkholderia
    pseudomallei* scored lowest among the three hits as indicated by its
    purplish hue, placing it near the bottom of the color-bar. However,
    the minimum value for the color-bar is approximately 80%. This
    result suggests that PanGIA’s confidence in correctly identifying
    this strain is somewhat lower than it’s confidence in correctly
    identifying the species or genus that the strain belongs to, but it
    still quite high. Dynamically scaled axes make comparison between
    close-scoring hits easier, given the landscape of hits are
    restricted to a limited range.

3.  The size of each hit on the scatterplot is scaled to its relative
    abundance in the sample. In this example, reads belonging to genus
    *Burkholderia* make up a comparatively much larger share than either
    species or strain.

4.  The associated data table for this PanGIA run contains three
    observations – corresponding to the three data points seen in the
    above scatterplot. A condensed summary of the information in a
    single row of this table appears as a tool-tip when the mouse is
    hovered over the center of a dot in the scatterplot.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/vis_datatable_example1.png">
<figcaption>
<em>the embedded datatable within the visualizer</em>
</figcaption>
</figure>

<br>

<div style="text-align: justify">

5.  The tool-tip presents numerical results from a given row of the data
    table. In the below image, *Burkholderia pseudomallei* is shown.

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/upclose_hovertool.png">
<figcaption>
<em>the tooltip appears when the mouse is over the center of a
datapoint</em>
</figcaption>
</figure>

<div style="text-align: justify">

## *Section 5)*: **Planned Visualizer Features**

The PanGIA GUI Visualizer is not complete. Two main elements – a
depth-scale-down chart and sample summary pie-charts – are still in
development. In addition, several widgets will be attached to the
scatterplot to make it more interactive, dynamic, and powerful. Planned
widget features include:

1.  Widgets for Numerical Filtering:
    1.  Read Count
    2.  Linear Coverage
    3.  Score
    4.  Depth-of-Coverage
    5.  Rank-Specific Metrics
2.  Widgets for User-Specified Filtering
    1.  Returning only specific genus/species/strain
    2.  Customize y-axis/size-axis/color-axis
3.  Widgets for Categorical Filtering:
    1.  Pathogenic Filtering
    2.  Rank Filtering

An example of the Visualizer as it will appear with widgets:

</div>

<figure>
<img src="gui/images_gui/Projects_and_Runs/scatterplot_with_widgets.png">
<figcaption>
<em>the final visualizer will include multiple sliders, dynamic axis
toggling, and other filtering tools</em>
</figcaption>
</figure>

## *Section 6)*: **Troubleshooting**

If having trouble with installation or usage, details in this section
may be of use.

### **Section 6.1)**: **Installation Troubleshooting**

We provide a second installation file, `gui_specs.txt` which may be of
use if the steps to construct a conda environment using a .yml file
specified in [**Section
3.1**](#section-31-preparing-a-non-dockerized-environment) do not work
correctly. Please follow these instructions:

1.  In a terminal window, run the command `conda deactivate` to ensure
    you begin in the base environment.
2.  Remove the environment you attempted to create in [**Section
    3.1**](#section-31-preparing-a-non-dockerized-environment) - if you
    used the default name ‘gui’, do this by running the command
    `conda env remove --name gui`. Skip this step if you’d prefer a
    different name for your working environment.
3.  In the same window, run the command
    `conda create --name gui python=3.6 --no-default-packages`. This
    creates a conda environment with none of the packages on the PATH of
    the base environment - ensuring that you have complete control of
    what is installed into your working environment.
4.  In the same window, run the command `conda activate gui` - replace
    ‘gui’ if you used a different name in Step 3.
5.  In the same window, run the command
    `conda install --name gui --file gui_specs.txt`.
6.  You should now be able to follow the instructions in [*Section
    4.1*](#section-411-non-dockerized-case) to launch the GUI.

### **Section 6.2)**: **eGPU Set-Up for GuPPY Basecalling**

If planning to use PanGIA as a Real-Time tool, using a GPU to power the
MinION/GuPPY software is required. This section covers a step-by-step
process for connecting an external GPU (eGPU) to a computer running the
Linux distribution of Pop-OS 20.04 - eventually, this section will be
expanded to cover set-up for Ubuntu 18.4.0 and 20.04 as well.

#### **Section 6.2.1)**: **Pop-OS 20.04**

These instructions were developed with a particular computer, OS, eGPU,
and eGFX chassis. As testing continues, we will consider other set-ups
as well. While using these instructions, mileage may vary for use with
other configurations. Please note however, that **the MinION/GuPPY
software is only compatible with Linux and plays most nicely with NVIDIA
GPUs**. We recommend downloading a NVidia GPU Driver to your computer
prior to following these instructions. The ingredients: 1. Computing
Unit: **NUC8:7HNK** 2. OS: **Pop-OS version 20.04** 3. eGPU: **NVidia
Quadro P4000** 1. NVidia GPU Driver: **460.73.01** 4. eGFX Chassis:
**EB3T-v3**

Connect the eGPU to the eGFX and ensure both are connected to their
power supply. Then, connect the eGFX to your computing unit via
Thunderbolt 3. Ensure that you are using a Thunderbolt cable and
connecting into Thunderbolt ports (look for the downward-arrow lightning
symbol). If the eGPU and eGFX chassis are compatible, the computer will
have access to the connected Thunderbolt device. You may
verify/troubleshoot this as follows: 1. In a terminal window, run the
command `cat /sys/bus/thunderbolt/devices/0-3/authorized`. If the
computer can see and has authorized access to the eGPU, this will return
a `1`. Please note that there may be variation in the above command
depending on your system - the `0-3` may instead be `0-0`, `0-1`, `0-2`,
or something else; the best way to verify is to navigate to
`/sys/bus/thunderbolt/devices/` and see what is present. 2. If the above
command does not return a `1`, run the command
`sudo sh -c 'echo 1 > /sys/bus/thunderbolt/devices/0-3/authorized'` in
order to grant access to the GPU, making any needed changes to the path
as described in Step 1. Verify eGPU visibility by running the command
`lspci | grep -i nvidia`, which will return information on the connected
eGPU.

Upgrade the Linux kernel and prepare your system for eGPU processing.
For each of these commands, press `y` when asked whether you’d like to
install/upgrade packages:

1.  `sudo apt-get update`
2.  `sudo apt-get upgrade`
3.  `sudo apt-get install gcc build-essential`
4.  `sudo apt-get install linux-headers-$(uname -r)`
5.  `sudo apt-get install linux-image-$(uname -r)`
6.  `sudo apt-get install linux-modules-extra-$(uname -r) linux-image-virtual`
7.  `reboot`

After the computer boots, it is time to install the NVidia driver -
**WARNING**: once the driver has been installed, there is a risk that
you will lose access to your OS’s graphical user interface. There are
detailed steps on how we fixed this problem using the above specified
set-up, but we can not guarantee it’s effectiveness for other cases.
That being said …

1.  If you haven’t already, download an NVidia GPU driver from the
    official site - **here’s the link** to the driver used in these
    instructions (**460.73.01**). We put ours in the Downloads folder.
2.  Blacklist the ‘nouveau’ graphics drivers that the eGPU will try to
    use by default. To do this:
    1.  Create the blacklist file by running command:
        `sudo touch /etc/modprobe.d/blacklist-nouveau.conf`.
    2.  Navigate to and open the newly created file via `nano` or some
        other command line editor.
    3.  Add two lines: `blacklist nouveau` and
        `options nouveau modeset=0` - save and close the file.

Next, update the config settings by running the following commands,
again pressing `y` when asked whether you’d like to update/upgrade
packages:

1.  `sudo update-initramfs -u`
2.  `sudo update-grub`
3.  `reboot`

Once the computer restarts, check whether the ‘nouveau’ drivers are
blacklisted by running the command `lsmod | grep -i nouveau` - this
should return nothing. From here, it is time to enter ‘Text-Only’ mode
to install NVidia drivers. Log out (without shutting down). Press
**Ctrl+Alt+F2** to enter ‘Text-Only’ mode, logging in with your username
and password. From here, you will disable the Linux service that powers
the OS’s GUI. Run the command `sudo service gdm3 stop`. If you are using
a different Linux OS, the name of the service may be different - instead
of `gdm3`, yours may be `lightdm`, `gdm`, or `sddm` - use whichever is
appropriate. You may verify that the service has been stopped by running
the command `service --status-all | grep -i` - look for a \[-\] next to
`gdm3`.

Drop a ‘run level’ by running the command `sudo init 3`. It may be a
good idea to check whether your eGPU is still recognized and authorized
by the computer, using aforementioned commands
`cat /sys/bus/thunderbolt/devices/0-3/authorized` and
`lspci | grep -i nvidia`.

Navigate to wherever you downloaded the NVidia driver and unzip/run it
using the command `sudo sh NIDIA-Linux-x86_64-460.73.01.run -x` -
replace the driver number (in this example, `460.73.01`) with whatever
driver you plan to use. Please note that the driver will be extracted
into the same folder.

Run the command `cd NVIDIA-Linux-x86_64-460.73.01/kernel/nvidia-uvm` to
navigate into the driver folder. Next:

1.  Navigate to and open the file named `uvm_va_block.c` with a command
    line editor. In your case, `uvm` may instead appear as `uvm8`.
2.  Add a line underneath the header composed of asterisks:
    `#include <linux/sched/task_stack.h>`.
3.  Look over the file - there should be a line of asterisks, followed
    by the line you just added, followed by the line
    `#include "uvm_linux.h`, which should have already been present.
4.  Save and close the file.

**Warning**: In this next phase, the driver installation will take
place. Please ensure that you have backed up everything you’d like to on
your hard drive - in the worst case, driver installation may make it
impossible to access your OS’s GUI. Proceed with caution. When you are
ready:

1.  Reboot your computer and press `F2` during start-up to open BIOS
    settings. In the Boot menu of the Advanced tab, find the box labeled
    **Secure Boot** and uncheck it. Boot the computer and re-enter Text
    Mode as described above (disable `gdm3` once more, if necessary).
2.  Navigate back to `/NVIDIA-Linux-x86_64-460.73.01` (as before,
    changing driver number as needed).
3.  Run the command: `sudo ./nvidia-installer --no-opengl-files` to
    initiate driver installation.
4.  You will be presented with a number of options during driver
    installation. We recommend the following parameters:
5.  Answer **Continue** or **Continue Anyway** when prompted about a
    pre-installation script failing.
6.  Answer **No** when prompted about **DKMS**.
7.  Answer **Yes** when asked about including 32-bit compatibility.
8.  Answer **No** when asked whether you’d like the new driver to
    **configure ‘X’**. This is especially crucial if you plan to use the
    eGPU only for computation, and never for accelerating the OS GUI’s
    graphics.
9.  Allow installation to complete.
