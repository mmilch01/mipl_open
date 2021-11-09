# mipl_open
Medical image processing and conversion tools developed at Neuroinformatics Lab Washington University St Louis

## Prerequisites
Prior to compilitng, download the <a href="https://dicom.offis.de">DCMTK 3.6.6</a> library and unzip to 'dcmtk' directory under the source root.

## Compiling
### Visual Studio 2019:
Use "mipl.sln" to open and build the entire solution.

### GCC 4.8.5 or higher:
edit 'mipl.def' to point SRCRT and RELEASE to sources and binary output dirs, respectively<br>

>make -f \<tool\>.mak [option]<br>
  
where \<tool\> is one of analyze2dcm, dcm_sort, dcminfo, imgcmp (former cmpanalyze).<br>
