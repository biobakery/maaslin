# MaAsLin: Multivariate Association with Linear Models

[MaAsLin](http://huttenhower.sph.harvard.edu/maaslin)

TBD

[MaAsLin Users Google Group](https://groups.google.com/forum/#!forum/maaslin-users)

MIT License

--------------------------------------------

## Contents ##
* [Description](#markdown-header-description)
* [Requirements](#markdown-header-requirements)
* [Installation](#markdown-header-installation)
* [Run a Demo](#markdown-header-run-a-demo)
* [How to Run](#markdown-header-how-to-run)
    * [Input Files](#markdown-header-input-files)
    * [Process Flow](#markdown-header-process-flow)
    * [Output Files](#markdown-header-output-files)
    * [Options](#markdown-header-options)
* [Troubleshooting](#markdown-header-troubleshooting)
* [How to Run in Galaxy](#markdown-header-how-to-run-in-galaxy)
* [Related Projects and Scripts](#markdown-header-related-projects-and-scripts)

## Description ##

MaAsLin is a multivariate statistical framework that finds
associations between clinical metadata and potentially
high-dimensional experimental data. MaAsLin performs boosted additive
general linear models between one group of data (metadata/the
predictors) and another group (in our case relative taxonomic
abundances/the response).  In our context we use it to discover
associations between clinical metadata and microbial community
relative abundance or function; however, it is applicable to other
data types.

Metagenomic data are sparse, and boosting is used to select metadata
that show some potential to be useful in a linear model between the
metadata and abundances. In the context of metadata and community
abundance, a sample's metadata is boosted for one Operational
Taxonomic Unit (OTU) (Yi). The metadata that are selected by boosting
are then used in a general linear model, with each combination of
metadata (as predictors) and OTU abundance (as response
variables). This occurs for every OTU and metadata combination. Given
we work with proportional data, the Yi (abundances) are
`arcsin(sqrt(Yi))` transformed. A final formula is as follows:

![](https://bitbucket.org/biobakery/maaslin/downloads/maaslinformula2.png)

## Requirements ##

MaAsLin requires the following R packages: agricolae, gam, gamlss, gbm, glmnet, inlinedocs, logging, MASS, nlme, optparse, outliers, penalized, pscl, robustbase

Please install these packages before installing MaAsLin.

## Installation ##

1. [Download](https://bitbucket.org/biobakery/maaslin/downloads/Maaslin_0.0.4.tar.gz) the latest version of MaAsLin.
2. Install MaAsLin (where X.Y.Z is the version number) : `` $ R CMD INSTALL Maaslin_X.Y.Z.tar.gz ``

## Run a Demo ##

Run the demo included in the MaAsLin install.

```
$ R
> library(Maaslin)
> example(Maaslin)
```

## How to Run ##

### Input Files ###
There are 2 essential input files: the "**\.read.config**" file and the "**input data**" file, and an optional "**\.R script**" . Details of each file follow:

** 1\. Input data File **

Required input file which we call the PCL file. This file that contains all the data
and metadata. This file is formatted so that metadata/data (otus or
bugs) are rows and samples are columns. All metadata rows should come
first before any abundance data. The file should be a tab delimited
text file. A demo PCL file is found in the MaAsLin download in maaslin/inst/extdata/ .

A PCL file is a TEXT delimited file similar to an excel spread sheet with the following characteristics. 

1. Rows represent metadata and features (bugs), columns represent samples.
2. The first row by default should be the sample ids.
3. Metadata rows should be next.
4. Lastly, rows containing features (bugs) measurements (like abundance) should be after metadata rows.
5. The first column should contain the ID describing the column. For metadata this may be, for example, "Age" for a row containing the age of the patients donating the samples. For measurements, this should be the feature name (bug name).
5. By default the file is expected to be TAB delimited.
6. If a consensus lineage or hierarchy of taxonomy is in the feature name, the default delimiter between clades is the pipe ("|").

** 2\. Read Config File **

Required input file \*.read.config. A read config file allows one to indicate what data is read from a PCL file without having to change the pcl file or change code. This means one can have a pcl file which is a superset of metadata and abundances which includes data you are not interested in for the run. This file is a text file with ".read.config" as an extension. This file is later described in detail in section **Process Flow ** subsection **4. Create your read.config file**.

** 3\. R Script File (Optional) **
 
Optional input file \*.R. The R script file is using a call back programming pattern that allows one to add/modify specific code to customize analysis without touching the main MaAsLin engine. A generic R script is provided "maaslin_demo2.R" and can be renamed and used for any study. The R script can be modified to add quality control or formatting of data, add ecological measurements, or other changes to the underlying data before MaAsLin runs on it. This file is not required to run MaAsLin.

### Updates in the new release ###
 

The program now can use as inputs either a PCL file or a TSV file.
It detects the type of the file by the suffix of the file  - *.tsv   or *.pcl

**Example of running with a tsv file: **

./R/Maaslin.R inst/extdata/maaslin_demo2.tsv output1  -i inst/extdata/maaslin_demo2.read.config

**Example of running with a pcl file: **

./R/Maaslin.R inst/extdata/maaslin_demo2.pcl output2 -i inst/extdata/maaslin_demo2.read.config


**The config file is optional - Can be generated dynamically **

If you pass the number of the row of the LAST Metadata (For a pcl file) or the number of the column of the LAST METADATA (For a tsv file) 
using the parameter --lastMetadata=  the Config file gets generated for you.

**Example of syntax to generate the config file automatically  for a pcl file:**

 ./R/Maaslin.R inst/extdata/maaslin_demo2.pcl output3  --lastMetadata=9
 
**Example of syntax to generate the config file automatically  for a tsv file:**
 
  ./R/Maaslin.R inst/extdata/maaslin_demo2.tsv output4  --lastMetadata=9



### Process Flow ###

** 1\. Obtain your abundance or relative function table. **

Abundance tables are normally derived from sequence data using
*Mothur*, *Qiime*, *HUMAnN*, or *MetaPhlAn*. Please refer to their documentation
for further details.

** 2\. Obtain your metadata. **

Metadata would be information about the samples in the study. For
instance, one may analyze a case / control study. In this study, you
may have a disease and healthy group (disease state), the sex of the
patents (patient demographics), medication use (chemical treatment),
smoking (patient lifestyle) or other types of data. All aforementioned
data would be study metadata. This section can have any type of data
(factor, ordered factor, continuous, integer, or logical
variables). If a particular data is missing for a sample for a
metadata please write NA. It is preferable to write NA so that, when
looking at the data, it is understood the metadata is missing and it's
absence is intentional and not a mistake. Often investigators are
interested in genetic measurements that may also be placed in the
metadata section to associate to bugs.

Please note that no special characters are allowed in the metadata
header names. These names should only contain alphanumeric 
characters with the addition of a period and an underscore: 
`[a-zA-Z0-9.-]`.

If you are not wanting to manually add metadata to your abundance
table, you may be interested in associated tools or scripts to help
combine your abundance table and metadata to create your pcl
file. Both require a specific format for your metadata file. Please
see the documentation for *QiimeToMaaslin* or *merge_metadata.py* (for
more details see section **Related Projects and Scripts**).

** 3\. Format and combine your abundance table and metadata as a pcl **

Please note two tools have been developed to help you! If you are
working from a Qiime OTU output and have a metadata text file try using
*QiimeToMaaslin* found at bitbucket. If you have a tab delimited file
which matches the below .pcl description (for instance MetaPhlAn
output) use the merge_metadata.py script provided in this project
(`maaslin/exec/merge_metadata.py`) and documented in
`maaslin/doc/Merge_Metadata_Read_Me.txt`.

** PCL format description **

i. Row 1 is expected to be sample IDs beginning the first column with a feature name to identify the row, for example "ID".

ii. Rows of metadata. Each row is one metadata, the first column entry
being the name of the metadata and each following column being the
metadata value for that each sample.

iii. Row of taxa/otu abundance. Each row is one taxa/otu, the first
column entry being the name of the taxa/otu followed by abundances of
the taxa/otu per sample.

iv. Abundances should be normalized by dividing each abundance measurement by the sum of the column (sample) abundances.  

v. Here is an example of the contents of an extremely small pcl file;
another example can be found in this project at
`maaslin/inst/extdata/maaslin_demo.pcl`.


    ID	Sample1	Sample2	Sample3	Sample4
    metadata1	True	True	False	False
    metadata2	1.23	2.34	3.22	3.44
    metadata3	Male	Female	Male	Female
    taxa1	    0.022	0.014	0.333	0.125
    taxa2	    0.406	0.029	0.166	0.300
    taxa3	    0.571	0.955	0.500	0.575


** 4\. Create your read.config file. **

A *.read.config file is a structured text file used to indicate which
data in the PCL file should be read into MaAsLin and used for
analysis. This allows one to keep their PCL file intact while
varying analysis. Hopefully, this avoids errors that may be introduced
while manipulating the PCL files.

Here is an example of the contents of a *.read.config file.

    Matrix: Metadata
    Read_PCL_Columns: Sample2-Sample15 
    Read_PCL_Rows: Age-Height,Weight,Sex,Cohort-Profession

    Matrix: Abundance
    Read_PCL_Columns: Sample2-Sample15
    Read_PCL_Rows: Bacteria-Bug100

The minimal requirement for a MaAsLin .read.config file is as
follows. The Matrix: should be specified. Metadata needs to be named
"Metadata" for the metadata section and "Abundance" for the abundance
section. “Read\_PCL\_Rows:” is used to indicate which rows are data or
metadata to be analyzed. Rows can be identified by their metadata/data
id. Separate ids by commas. If there is a consecutive group of
metadata/data a range of rows can be defined by indicating the first
and last id separated by a “-“. If the beginning or ending id is
missing surrounding an “–“, the rows are read from the beginning or to
the end of the pcl file, respectively.

A minimal example is shown here:

    Matrix: Metadata
    Read_PCL_Rows: -Weight

    Matrix: Abundance
    Read_PCL_Rows: Bacteria-

With this minimal example, the delimiter of the file is assumed to be
a tab, all columns are read (since they are not indicated
here). Metadata are read as all rows from the beginning of the pcl
file (skipping the first Sample ID row) to Weight; all data are read
as all rows from Bacteria to the end of the pcl file. This example
refers to the default input files given in the MaAsLin download as
maaslin/inst/extdata/maaslin_demo2.\*.

** 5\. Create your R script (Optional). **

The R script is used to add code that manipulates your data before
analysis, and for manipulating the multifactoral analysis figure. A
default “*.R” script is available with the default MaAsLin project at
maaslin/inst/extdata/maaslin_demo2.R. This is an expert option and should
only be used by someone very comfortable with the R language.

** 6\. Specify Output Directory **

If you have a specific directory where the results must save, then specify it's location. If not, specify the name of the directory to create with it's path. E.g. If you want results to go the folder "Results" in your current path, specify "Results". If it must go to the directory "Results" in another location, specify "Path of other location/Results".

** 7\. Run. **

If starting with a PCL file (input.pcl), first transpose it to a TSV file (input.tsv).

`` $ ./Maaslin/exec/transpose.py < input.pcl > input.tsv ``

Run MaAsLin.

`` $ ./R/Maaslin.R -i input.read.config input.tsv outputdir ``
 ```
```
$ R
> library(Maaslin)
> Maaslin('input.tsv','maaslin_output',strInputConfig='input.read.config')
```

Please see the [FAQs](https://bitbucket.org/biobakery/maaslin/wiki/Home) if you need information on running MaAsLin from the command line.

** 8\. Discover amazing associations in your results! **


### Output Files ###

The following files will be generated per MaAsLin run. In the
following listing the term $PROJECTNAME refers to what you named your PCL file without the extension. There are four main types of output files. Additional files are created in a folder named "QC".

**1\. $PROJECTNAME_log.txt **

This file contains the detail for the statistical engine. This can be
useful for detailed troubleshooting.

**2\. $PROJECTNAME-$METADATA.txt**

Each metadata will have a file of associations. $METADATA
is the metadata column name in the PCL file. Any associations
indicated to be performed after initial variable selection (boosting)
is recorded here. Included are the information from the final general
linear model (performed after the boosting) and the FDR corrected
p-value (q-value). Can be opened as a text file or spreadsheet.

**3\. $PROJECTNAME-$METADATA.pdf**

Any association that had a q-value less than or equal to the given
significance threshold will be plotted here (default is 0.25; can be
changed using the commandline argument -d). If this file does not
exist, the $PROJECTNAME-$METADATA.txt should not have an entry that is
less than or equal to the threshold. Factor data is plotted as
knotched box plots; continuous data is plotted as a scatter plot with
a line of best fit. 

**4\. $PROJECTNAME.txt**

A collection of all entries in the $PROJECTNAME-$METADATA.pdf. Can be
opened as a text file or spreadsheet.

**5\. QC/data.tsv **

The data matrix that was read in (transposed). Useful for making sure
the correct data was read in.

**6\. QC/data.read.config **

Can be used to read in the data.tsv.

**7\. QC/metadata.tsv **

The metadata that was read in (transposed). Useful for making sure the
correct metadata was read in.

**8\. QC/metadata.read.config **

Can be used to read in the data.tsv.

**9\. QC/read-Merged.tsv **

The data and metadata merged (transposed). Useful for making sure the
merging occurred correctly.

**10\. QC/read-Merged.read.config **

Can be used to read in the read_merged.tsv.

**11\. QC/read_cleaned.tsv **

The data read in, merged, and then cleaned. After this process the
data is written to this file for reference if needed.

**12\. QC/read_cleaned.read.config **

Can be used to read in read_cleaned.tsv.

**13\. QC/ProcessQC.txt **

Contains quality control for the MaAsLin analysis. This includes
information on the magnitude of outlier removal.

**14\. QC/Run_Parameters.txt **

Contains an account of all the options used when running MaAsLin so the exact methodology can be recreated if needed.

### Options ###

Although we recommend the use of default options, arguments exist to modify both MaAsLin methodology and figures. To see
an up-to-date listing of argument usage, type ` > help(Maaslin)`.

**Example args:**

1. strVerbosity = "DEBUG" 
1. dSignificanceLevel = 0.1 

In this example MaAsLin is modified to produce verbose output for
debugging and to change the threshold for making pdfs to a
q-value equal to or less than 0.1.

** All verses All **

The all verses all analysis flow is a way of manipulating how metadata are used. In this method there is a group of metadata that are always evaluated, as well there are a group that are added to this one at a time. To give a more concrete example: You may have metadata cage, diet, and treatment. You may always want to have the association of abundance evaluated controlling for cage but otherwise looking at the metadata one at a time. In this way the cage metadata is the forced part of the evaluation while the others are not forced and evaluated in serial. The appropriate options to indicate this follows:

1. fAllvAll = TRUE
1. strForcedPredictors = "cage"

(fAllvAll=TRUE) indicates all verses all is being used. 

** Forcing covariates **
(strForcedPredictors) indicates which metadata are forced (multiple metadata can be given comma delimited as shown here =c("metadata1","metadata2","metadata3")). This does not bypass the feature selection method so the metadata that are not forced are subject to feature selection and may be removed before coming to the evaluation. If you want all the metadata that are not forced to be evaluated in serial you will need to turn off feature selection and will have final combined options as seen here:

1. fAllvAll = TRUE
1. strForcedPredictors = "cage"
1. strModelSelection = "none"

## Troubleshooting ##

Please see the [FAQs](https://bitbucket.org/biobakery/maaslin/wiki/Home) for MaAsLin troubleshooting.

## How to Run in Galaxy ##

MaAsLin can be run in galaxy at : http://huttenhower.sph.harvard.edu/galaxy/

Please refer to the Huttenhower Galaxy site for details on how to run MaAsLin in galaxy.

## Related Projects and Scripts ##
Other projects exist that may help in your analysis:

** QiimeToMaAsLin **
    QiimeToMaAsLin is a project that reformats abundance files from
    Qiime for MaAsLin. Several formats of Qiime consensus lineages are
    supported for this project. To download please visit
    [https://bitbucket.org/biobakery/qiimetomaaslin](https://bitbucket.org/biobakery/qiimetomaaslin).

** Merge_Metadata ** 
    merge_metadata.py is a script included in the MaAsLin project to
    generically merge a metadata file with a table of microbial (or
    other) measurements. This script is located in `maaslin/exec` and
    is documented in `maaslin/doc/ Merge_Metadata_Read_Me.txt`
