There are 2 files that read the RES and SPX files and store run information:

1) MfixData.cpp
2) MfixData.h


Included is a test routine:

mfix_reader.cpp


To test the code out ...

1) edit the test routine, mfix_reader.cpp

   Near the beginning of main() is a statement:

      data.SetName("BUB01");

   This would be for BUB01.RES , BUB01.SP1 , etc.

   Change the "BUB01" to your run name


2) compile/link .. example using g++

   g++ mfix_reader.cpp MfixData.cpp -o mfix_reader


3) run the code

   ./mfix_reader


4) Look at the output in run_info.txt




a) General usage ... required steps


   #include "MfixData.h"

   // ...

   MfixData data;

   data.SetName("BUB01");  // change as needed
     
   data.ReadRes0();
   
   data.CreateVariableNames();

   data.GetTimes();        // loops thru all SPx files, finds
                           // what times are in the files and
                           // the position in the files that
                           // they are at



b) Look at the run_info.txt output file and the corresponding code.
   That should give you a start.  You can contact me with any
   questions that come up.

c) Some examples:

   1) dimensions in the variables : data.imax2 , data.jmax2 ,
                                    data.kmax2 , data.ijkmax2


   2) Currently, the code does not have arrays for all the field variables
      (I was using the code to compare two different large runs, so I kept
      memory to a minimum).

   
   3) The following loop would iterate over all the times and read EP_g.

   int variable_index = 0;       // EP_g ... see run_info.txt file for
                                 //          a mapping between index and variable

   int ijk = 234;  // just some cell for sample output

   for (int ts = 0; ts<data.Get_NTIMES()-1; ++ts)
   {
      data.GetVariableAtTimestep(var_index,ts);

      // reads into std::vector<float> data.var

      out << "The value of " << data.variable_names[var_index]
          << "[" << ijk << "] = " << data.var[ijk]
          << " at timestep index = " << ts << "\n";
   }


   4) So currently, if you wanted the values for both EP_g and P_g you would
      need to do something like this:

      data.GetVariableAtTimestep(0,ts);
      std::vector<float> EP_g = var;

      data.GetVariableAtTimestep(1,ts);
      std::vector<float> P_g = var;


   5) To change it so that all the field variables are declared so that
      the above copies are not needed:

      a) you would need to declare the variables in MfixData class (MfixData.h).

      b) Change the GetVariableAtTimestep member function (it is the last
         function in MfixData.cpp)

      c) If you want, I can help you with this.


      I can help you if you think the code can be modified to fit your needs.


D) If you look at the sample_run_info.txt file, you will notice that EP_g
   was written to the SPx file every 0.01 seconds, while P_g was only written
   every 0.1 seconds.  If you request P_g at a time = 0.55 , it returns the
   value at time = 0.5 (it does not do any interpolation).
