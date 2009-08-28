#include "MfixData.h"

#include <fstream>  // for debugging pursposes
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

int main()
{
    int i , ijk;

    ofstream out("run_info.txt");  // for debugging purposes

    MfixData data;

    data.SetName("BUB01");
     
    data.ReadRes0();
   
    data.CreateVariableNames();

    data.GetTimes();
   

    // rest is for debugging / information
   
    out << " size : " ;
    out << data.imax2 << " * " << data.jmax2 << " * "
        << data.kmax2 << " = " << data.ijkmax2 << "\n";

	//////////////////////////////////////////////////////////////////////////
    //                        output variable names                         //

    out << "There are " << data.Get_NVARS() << " variables defined\n\n\t";
   
    for (i=0; i<data.variable_names.size(); ++i)
    {
       out << data.variable_names[i] << "\n\t";
    }

	//////////////////////////////////////////////////////////////////////////
    //                        output time step data

    out << "\n\n\nThere are " << data.Get_NTIMES() << " times\n\n\t";


    for (i=0; i<data.times.size(); ++i)
    {
       out << data.times[i] << "\n\t";
    }

	//////////////////////////////////////////////////////////////////////////
    //           output timestep / file position information                //
   
    out << "\n\n";

    char ext [] = "123456789ab";

    for (i=0; i<data.nspx_use; ++i)
    {
        out << data.run_name << ".SP" << ext[i] << "\n";

        map<float,ULONGLONG> & tmap = data.timeMap2[i];
        map<float,ULONGLONG>::iterator it;
        for (it=tmap.begin(); it!=tmap.end(); ++it)
        {
            out << "\ttime = " << it->first << " at position = " << it->second << "\n";
        }

    }

	//////////////////////////////////////////////////////////////////////////
    //             output ... paraview index variable map to MfixData       //
	//                        internal variable index                       //

    out << "\n\n" << setw(10) << "Paraview" 
                  << setw(10) << "MfixData"  << "\t"
                  << setw(10) << "variable*"
                  << "\n\n";

    for (int j=0; j<data.index_map.size(); ++j)
    {
        out << setw(10) << j 
            << setw(10) << data.index_map[j] << "\t"
            << setw(10) << data.variable_names[ data.index_map[j] ]
            << "\n";
    }

	//////////////////////////////////////////////////////////////////////////
    //                 output value of variable at a time step              //

    int var_index = 0; //MfixData index , not paraview index

    ijk = 15; 
   
    for (int ts = 0; ts<data.Get_NTIMES(); ++ts)
    {

       data.GetVariableAtTimestep(var_index,ts);

       out << "The value of " << data.variable_names[var_index]
           << "[" << ijk << "] = " << data.var[ijk]
           << " at timestep index = " << ts << "\n";
    }

	//////////////////////////////////////////////////////////////////////////
   
    out << "\n\n       Flags                    \n\n";
   
    vector<int> vPOUT;
   
    const int flag_lookup = 20;
    map<int,int> flag_freq;
   
    for (ijk=0; ijk<data.ijkmax2; ++ijk)
    {
       ++flag_freq[ data.FLAG[ijk] ];
       if (data.FLAG[ijk] == flag_lookup) vPOUT.push_back(ijk);
    }
   
    //
   
    out << " flag frequency ... \n";
    map<int,int>::iterator it = flag_freq.begin();
    while (it != flag_freq.end())
    {
       out << "\t" << it->first << "\t:" << it->second << "\n";
       ++it;
    }
     
    out << "\n\ncells with flag = " << flag_lookup << "\n";
    for (ijk=0; ijk<vPOUT.size(); ++ijk)
    {
       out << "\t" << vPOUT[ijk] << "\n";
    }
   
	//////////////////////////////////////////////////////////////////////////


    return 0;
}
