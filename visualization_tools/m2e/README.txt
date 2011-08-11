
MFIX to Ensight Translator : m2e

1) to create the translator

   make

2) to convert MFIX RES and SPx files to Ensight

   a. copy the m2e executable to the folder with the MFIX output files
   b. ./m2e RUNNAME

      (i.e. enter the name of the RES file without the .RES)

3) Bring up Ensight

   a. Choose the "File" ... "open" menu

   b. In the dialog,

      1. change the "format" to "Case"
      2. Highlight the RUNNAME.cas file
      3. click on the "set case" button
      4. click on the "load all" button

   c. The part you would normally choose to plot is:

      "grid for cell data"



NOTE: On later versions of Ensight, some openGL features might
      not work when running Ensight remotely.  If everything
      seems to read correctly, but no image is displayed,
      you can try setting the following before running
      Ensight ....

      setenv CEI_NUM_AUX_BUFFERS 0

