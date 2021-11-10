# mipl_open
Medical image processing and conversion tools developed at Neuroinformatics Lab Washington University St Louis

## Prerequisites
Prior to compilitng, download and build the <a href="https://dicom.offis.de">DCMTK 3.6.6</a> library. Put the build in 'dcmtk' directory under the source root.

## Compiling
### Visual Studio 2019:
Use "mipl.sln" to open and build the entire solution.

### GNU Make v. 3.82+, GCC 4.8.5+:
>make -f \<tool\>.mak [clean|cleanall]<br>
  
where \<tool\> is one of analyze2dcm, dcm_sort, dcminfo, imgcmp (former cmpanalyze).<br>

