#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cstring>
#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

#include "MfixData.h"

// MFIX restart helper functions


namespace
{
    void SWAP_DOUBLE(double & value)
    {
    
        static char Swapped[8];

        double * Addr = &value;

        Swapped[0]=*((char*)Addr+7);
        Swapped[1]=*((char*)Addr+6);
        Swapped[2]=*((char*)Addr+5);
        Swapped[3]=*((char*)Addr+4);
        Swapped[4]=*((char*)Addr+3);
        Swapped[5]=*((char*)Addr+2);
        Swapped[6]=*((char*)Addr+1);
        Swapped[7]=*((char*)Addr  );

        value = *(reinterpret_cast<double*>(Swapped));
    
    }

    void SWAP_FLOAT(float & value)
    {
    
        static char Swapped[4];

        float * Addr = &value;

        Swapped[0]=*((char*)Addr+3);
        Swapped[1]=*((char*)Addr+2);
        Swapped[2]=*((char*)Addr+1);
        Swapped[3]=*((char*)Addr  );

        value = *(reinterpret_cast<float*>(Swapped));
    
    }

    void SWAP_INT(int & value)
    {
    
        static char Swapped[4];

        int * Addr = &value;

        Swapped[0]=*((char*)Addr+3);
        Swapped[1]=*((char*)Addr+2);
        Swapped[2]=*((char*)Addr+1);
        Swapped[3]=*((char*)Addr  );

        value = *(reinterpret_cast<int*>(Swapped));
    
    }


    string ConvertIntToString(int in)
    {
        string out;
        stringstream buff;
        buff << in;
        buff >> out;
        return out;
    }

    int ConvertCharToInt(char in)
    {
        int out;
        stringstream buff;
        buff << in;
        buff >> out;
        return out;
    }

    int ConvertStringToInt(const string & in)
    {
        int out;
        stringstream buff;
        buff << in;
        buff >> out;
        return out;
    }


    void GetInt(istream& in, int &val)
    {
        in.read( (char*)&val,sizeof(int));
        SWAP_INT(val);
    }

    void GetDouble(istream& in, double& val)
    {
        in.read( (char*)&val,sizeof(double));
        SWAP_DOUBLE(val);
    }

    void GetFloat(istream& in, float& val)
    {
        in.read( (char*)&val,sizeof(float));
        SWAP_FLOAT(val);
    }

    static char ext [] = "123456789AB";

    char buffer[513];

    std::string version;
    double version_number;

    double P_ref;
    double P_scale;
    int    DIM_IC;
    int    DIM_BC;
    int    DIM_C;
    int    DIM_IS;
    int    next_reca;

    double c_e;
    double c_f;
    double phi;
    double phi_w;

    double C_e, C_f, Phi, Phi_w;

    int spx_records_per_timestep;

    map<string,int> varNameToSPX;
    vector<double> C;
   
    double dt , xmin ;

    void SkipBytes(istream& in, int n)
    {
       in.read(buffer,n); // maybe seekg instead
    }
}


  

MfixData::MfixData()
{
   P_ref     = 0;
   P_scale   = 1;
   DIM_IC    = 5;
   DIM_BC    = 5;
   DIM_C     = 5;
   DIM_IS    = 5;
   c_e       = 1.0;
   c_f       = 0.0;
   phi       = 0.0;
   phi_w     = 0.0;
   nspx_use  = 9;
   NScalar   = 0;
   bKepsilon = false;  
   nRR       = 0;
   nTimes    = 0;

   res_option = 0;
}

MfixData::~MfixData() { }



void MfixData::RestartVersionNumber(char* buffer)
{
   string s;

   stringstream ss(buffer);
   ss >> s >> s >> version_number;  // RES = 01.10
   version = buffer;
   version.resize(11);
}


void MfixData::IN_BIN_512(istream& in, double* v, int n)
{
   const int nr = 512/sizeof(double);

   double array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( reinterpret_cast<char*>(&array) , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              v[c] = array[j];
              SWAP_DOUBLE(v[c]);
              ++c;
           }
       }
   }
}

void MfixData::IN_BIN_512R(istream& in, float* v, int n)
{
   const int nr = 512/sizeof(float);
   float array[nr];

   int c = 0;
   for (int i=0; i<spx_records_per_timestep; ++i)
   {
       in.read( reinterpret_cast<char*>(&array) , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              v[c] = array[j];
              SWAP_FLOAT(v[c]);
              ++c;
           }
       }
   }
}

void MfixData::IN_BIN_512I(istream& in, int* v, int n)
{
   const int nr = 512/sizeof(int);

   int array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( reinterpret_cast<char*>(&array) , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              v[c] = array[j];
              SWAP_INT(v[c]);
              ++c;
           }
       }
   }
}

void MfixData::ReadRes0()
{
   variable_names.clear();
   varNameToSPX.clear();
   index_map.clear();
   timeMap2.clear();
   times.clear();

   int i,l,n,lc;
    
   // init some variables that are not read in in older versions
   
   int DIMENSION_USR = 5;
   
   string fname = run_name + ".RES";
//   fstream in(fname.c_str(),ios::binary|std::ios::in|std::ios::out);
   ifstream in(fname.c_str(),ios::binary|std::ios::in);

   if (!in)
   {
      cout << "could not open file" << endl;
      return;
   }

   // version : record 1
   memset(buffer,0,513);
   in.read(buffer,512);
   RestartVersionNumber(buffer);
   

   // skip 2 lines : records 2 and 3

   if (res_option != 3626262)  // was == 3
   {
       in.read(buffer,512);
   }
   else
   {
      /*
      size_t pos_rec2 = in.tellg();
      in.read(buffer,60);

      int id_month , id_day , id_year , id_hour , id_minute , id_second;
      in.read( reinterpret_cast<char*>(&id_month),sizeof(int));
      in.read( reinterpret_cast<char*>(&id_day),sizeof(int));
      in.read( reinterpret_cast<char*>(&id_year),sizeof(int));
      in.read( reinterpret_cast<char*>(&id_hour),sizeof(int));
      in.read( reinterpret_cast<char*>(&id_minute),sizeof(int));
      in.read( reinterpret_cast<char*>(&id_second),sizeof(int));

      SWAP_INT(id_month);
      SWAP_INT(id_day);
      SWAP_INT(id_year);
      SWAP_INT(id_hour);
      SWAP_INT(id_minute);
      SWAP_INT(id_second);

      cout << "month  in RES file = " << id_month << "   enter new value > ";
      cin  >> id_month;

      cout << "day    in RES file = " << id_day   << "   enter new value > ";
      cin  >> id_day;

      cout << "year   in RES file = " << id_year  << "   enter new value > ";
      cin  >> id_year;

   
      cout << "hour   in RES file = " << id_hour  << "   enter new value > ";
      cin  >> id_hour;

      cout << "minute in RES file = " << id_minute << "   enter new value > ";
      cin  >> id_minute;

      cout << "second in RES file = " << id_second << "   enter new value > ";
      cin  >> id_second;
   
   

      SWAP_INT(id_month);
      SWAP_INT(id_day);
      SWAP_INT(id_year);
      SWAP_INT(id_hour);
      SWAP_INT(id_minute);
      SWAP_INT(id_second);

      in.seekp(pos_rec2);
      in.write(buffer,60);

      in.write( reinterpret_cast<char*>(&id_month),sizeof(int));
      in.write( reinterpret_cast<char*>(&id_day),sizeof(int));
      in.write( reinterpret_cast<char*>(&id_year),sizeof(int));
      in.write( reinterpret_cast<char*>(&id_hour),sizeof(int));
      in.write( reinterpret_cast<char*>(&id_minute),sizeof(int));
      in.write( reinterpret_cast<char*>(&id_second),sizeof(int));

      cout << "\ndates changed\n\n";

      return;
      */
   }




   in.read(buffer,512);

   // imin1 etc : record 4
   memset(buffer,0,513);
   
   if (version == "RES = 01.00")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetDouble(in,dt);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 15 ints ... 4 doubles = 92 bytes
      SkipBytes(in,420);
   }
   else if(version == "RES = 01.01" || version == "RES = 01.02")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC);    GetInt(in,DIM_BC);
      GetDouble(in,dt);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 17 ints ... 4 doubles = 100 bytes
      SkipBytes(in,412);
   }
   else if(version == "RES = 01.03")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 17 ints ... 5 doubles = 108 bytes
      SkipBytes(in,404);
   }
   else if(version == "RES = 01.04")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 18 ints ... 5 doubles = 112 bytes
      SkipBytes(in,400);
   }
   else if(version == "RES = 01.05")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetInt(in,DIM_IS);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 19 ints ... 5 doubles = 116 bytes
      SkipBytes(in,396);
   }
   else
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetInt(in,DIM_IS);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
      GetDouble(in,C_e); GetDouble(in,C_f); GetDouble(in,Phi); GetDouble(in,Phi_w);
      // 19 ints ... 9 doubles = 148 bytes
      SkipBytes(in,364);
   }

   cout << "imax2 = " << imax2 << "\n";
   cout << "jmax2 = " << jmax2 << "\n";
   cout << "kmax2 = " << kmax2 << "\n";
   
   const int nr = 512/sizeof(float);

   if ( ijkmax2%nr == 0)
      spx_records_per_timestep = ijkmax2/nr;
   else
      spx_records_per_timestep = 1 + ijkmax2/nr;
      
 
   var.resize(ijkmax2);

   // C , C_name and nmax

   NMAX.resize(MMAX+1);
   for (lc=0; lc<MMAX+1; ++lc) NMAX[lc] = 1;
   
   C.resize(DIM_C);

   if (version_number > 1.04)
   {
      IN_BIN_512 (in,&C[0],DIM_C);

      for (lc=0; lc<DIM_C; ++lc) in.read(buffer,512);  // c_name[] 
 
      if (version_number < 1.12)
          IN_BIN_512I(in,&NMAX[0],MMAX+1);
      else
      {  
          size_t pos = in.tellg();

          IN_BIN_512I(in,&NMAX[0],MMAX+1);

          if (res_option == 11111111)  // was == 1
          {
             /*
             cout << "\n\nif imax2/jmax2/kmax2 are incorrect ... kill program\n\n";

             cout << "MMAX = " << MMAX << "\n";

             for (int qq=0; qq<MMAX+1; ++qq)
             {
                 cout << "value of NMAX(" << qq << ") in RES file = " << NMAX[qq]
                      << "    ...  enter new value > ";

                 cin >> NMAX[qq];

                 SWAP_INT( NMAX[qq] );
             }
             cout << "attempting ...\n";
             if (!in) cout << "problem at start\n";
             in.seekp(pos);
             in.write( reinterpret_cast<char *>( &NMAX[0] ), (MMAX+1)*sizeof(NMAX[0]) );
             if (!in) 
                cout << "problem at end\n";
             else
                cout << "file has been modified\n";

            return;
            */
         }
      }
   }
  
   DX.resize(imax2);
   DY.resize(jmax2);
   DZ.resize(kmax2);

   IN_BIN_512(in,&DX[0],imax2);
   IN_BIN_512(in,&DY[0],jmax2);
   IN_BIN_512(in,&DZ[0],kmax2);
   
 
   // run_name etc.
   
   memset(units,0,17);
   memset(coordinates,0,17);
   
   in.read(buffer,120);      // run_name , description
   in.read(units,16);        // units
   in.read(buffer,16);       // run_type
   in.read(coordinates,16);  // coordinates 
   SkipBytes(in,512-168);

   char tmp[17];
   
   memset(tmp,0,17);

   int ic = 0;
   for (i=0; i<17; ++i)
   {
        if (units[i] != ' ') tmp[ic++] = units[i];
   }

   memset(tmp,0,17);

   ic = 0;
   for (i=0; i<17; ++i)
   {
        if (coordinates[i] != ' ') tmp[ic++] = coordinates[i];
   }
   strcpy(coordinates,tmp);

   vector<int>    tmpI;
   vector<double> tmpD;

   if (version_number >= 1.04)
   {
      tmpD.resize(NMAX[0]);
      IN_BIN_512(in,&tmpD[0],NMAX[0]);             // MW_g
      for (i=0; i<MMAX; ++i) in.read(buffer,512);  // MW_s
   }
   in.read(buffer,512);  // D_p etc.

   // read in the "DIM_IC" variables (and ignore ... not used by ani_mfix)
   tmpI.resize(DIM_IC);
   tmpD.resize(DIM_IC);

   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_x_w
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_x_e
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_y_s
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_y_n
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_z_b
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_z_t ...
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_i_w
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_i_e
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_j_s
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_j_n
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_k_b
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_k_t ...
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_ep_g
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_p_g
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_t_g
   
//   cout << "after ic_t_g " << in.tellg() << "\n";
//   cout << "vn = " << version_number << "\n";

   if (version_number < 1.15)
   {
      IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_t_s(1,1)
      IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_t_s(1,2) or ic_tmp 
   }

   if (version_number >= 1.04)
   {
      for (int i=0; i<NMAX[0]; ++i) IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_x_g
   }
//   cout << "after ic_x_g " << in.tellg() << "\n";

   IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_u_g
   IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_v_g
   IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_w_g

   for (lc=0; lc<MMAX; ++lc)
   {
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_rop_s
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_u_s
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_v_s
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_w_s
      
 //     cout << "after ic_w_s " << in.tellg() << "\n";
      if (version_number >= 1.15)
      {
         IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_t_s
      }
      
      if (version_number >= 1.04)
      {
         for (n=0; n<NMAX[lc+1]; ++n) // ...
            IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_x_s
      }
   }

   // read in the "DIM_BC" variables (and ignore ... not used by ani_mfix)
   tmpI.resize(DIM_BC);
   tmpD.resize(DIM_BC);
   
//   cout << "bc_start = " << in.tellg() << "\n";


   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_x_w
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_x_e
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc y s
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc y n
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc z b
   IN_BIN_512(in,&tmpD[0],DIM_BC);  // bc z t
   IN_BIN_512I(in,&tmpI[0],DIM_BC);  // bc i w
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc i e
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc j s
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc j n
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc k b
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc k t
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc ep g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc p g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc t g

   if (version_number < 1.15)
   {
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_t_s(1,1)
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_t_s(1,1) or bc_tmp
   }

   if (version_number >= 1.04)
   {
      for (int i=0; i<NMAX[0]; ++i) IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_x_g
   }

   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc u g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc v g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc w g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc ro g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_rop_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc volflow g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc massflow g

   for (lc=0; lc<MMAX; ++lc)
   {
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc rop s
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc u s
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc v s
      
      if (version_number >= 1.04)
      {
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc w s

         if (version_number >= 1.15)
         {
            IN_BIN_512(in,&tmpD[0],DIM_BC); // bc t s
         }
         for (n=0; n<NMAX[lc+1]; ++n) // ...
         {      
            IN_BIN_512(in,&tmpD[0],DIM_BC); // bc x s
         }
      }
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc volflow s
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc massflow s
   }


   if (version == "RES = 01.00")
      l = 10;
   else
      l = DIM_BC;

   for (lc=0; lc<l; ++lc) in.read(buffer,512); // BC TYPE
   
   size_t flag_pos = in.tellg();
//   cout << " before flag " << in.tellg() << "\n";

   FLAG.resize(ijkmax2);
   IN_BIN_512I(in,&FLAG[0],ijkmax2);

   if (res_option == 2222222)  // was == 2
   {
      /*
      int cccc = 0;
      for (int qqq=0; qqq<ijkmax2; ++qqq) 
      {
         if (FLAG[qqq] == 1) cccc = cccc + 1;
         FLAG[qqq] = 1;
         SWAP_INT( FLAG[qqq] );
      }
      cout << "number of initial FLAG=1 cells = " << cccc << " out of " << ijkmax2 << " cells \n";

      in.seekp(flag_pos);
      in.write( reinterpret_cast<char*>(&FLAG[0]),ijkmax2*sizeof(FLAG[0]));

      cout << "\nfile modified\n";

      return;
      */
   }
      

   // DIM_IS varibles (not needed by ani_mfix)
   tmpI.resize(DIM_IS);
   tmpD.resize(DIM_IS);

   if (version_number >= 1.04)
   {
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is x w
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is x e
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is y s
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is y n
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is z b
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is z t
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is i w
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is i e
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is j s
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is j n
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is k b
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is k t
      IN_BIN_512(in,&tmpD[0],DIM_IS);  // is_pc(1,1)
      IN_BIN_512(in,&tmpD[0],DIM_IS);  // is_pc(1,2)
     
      if (version_number >= 1.07)
      {
         for (l=0; l<MMAX; ++l) IN_BIN_512(in,&tmpD[0],DIM_IS); // is_vel_s
      }

      for (lc=0; lc<DIM_IS; ++lc) in.read(buffer,512); // is_type
   }

   if (version_number >= 1.08) in.read(buffer,512);
   
//   cout << "after >= 1.08 " << in.tellg() << "\n";
//   
   if (version_number >= 1.09) 
   {
      in.read(buffer,512);
      
      if (version_number >= 1.5)
      {
         GetInt(in,nspx_use);
         SkipBytes(in,508);
      }
 //     cout << "after n_spx_use " << in.tellg() << "\n";
 //     cout << nspx_use << "\n";
      
      for (lc=0; lc< nspx_use; ++lc) in.read(buffer,512); // spx_dt
      
      for (lc=0; lc<MMAX+1; ++lc) in.read(buffer,512);    // species_eq
      
      tmpD.resize(DIMENSION_USR);
      
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr x w
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr x e
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr y s
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr y n
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr z b
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr z t
     
 //     cout << "after usr z t " << in.tellg() << "\n";
      for (lc=0; lc<DIMENSION_USR; ++lc) in.read(buffer,512);    // usr_ext etc.
        
          
      tmpD.resize(DIM_IC);      
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_p_star
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_l_scale
      for (lc=0; lc<DIM_IC; ++lc) in.read(buffer,512);    // ic_type
          
      tmpD.resize(DIM_BC);      
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_dt_0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_jet_g0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_dt_h
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_jet_gh
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_dt_l
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_jet_gl
 //     cout << "after bc_jet_gl " << in.tellg() << "\n";
   }

      
   if (version_number >= 1.1)  in.read(buffer,512);  // mu_gmax
   if (version_number >= 1.11) in.read(buffer,512);  // x_ex , model_b
      
   if (version_number >= 1.12)
   {
      in.read(buffer,512);   // p_ref , etc.
      in.read(buffer,512);   // leq_it , leq_method
    
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_hw_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_uw_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_vw_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_ww_g
    
      for (lc=0; lc<MMAX; ++lc)
      {
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_hw_s
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_uw_s
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_vw_s
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_ww_s
      }
   }
//   cout << "after bc_ww_s " << in.tellg() << "\n";
      
   if (version_number >= 1.13) in.read(buffer,512);    // momentum_x_eq , etc.
   if (version_number >= 1.14) in.read(buffer,512);    // detect_small
      
   if (version_number >= 1.15)
   {
      in.read(buffer,512);    // k_g0 , etc.
     
      tmpD.resize(DIM_IC);   
     
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_gama_rg
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_t_rg
        
      for (lc=0; lc<MMAX; ++lc)
      {
         IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_gama_rs
         IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_t_rs
      }
  }
     
  if (version_number >= 1.2) in.read(buffer,512); // norm_g , norm_s
 
  if (version_number >= 1.3)
  {
     GetInt(in,NScalar);
     SkipBytes(in,sizeof(double)); // tol_resid_scalar

     int DIM_tmp;
     GetInt(in,DIM_tmp);
     SkipBytes(in,512-sizeof(double)-2*sizeof(int));
    
     tmpI.resize(DIM_tmp);
     IN_BIN_512I(in,&tmpI[0],DIM_tmp);  // Phase4Scalar;
  }
        
  if (version_number >= 1.5)
  {
     GetInt(in,nRR);
     SkipBytes(in,508);
  }
     
     
  if (version_number >= 1.5999)
  {
     int tmp;
     GetInt(in,tmp);
     SkipBytes(in,508);
    
     if (tmp != 0) bKepsilon = true;
  }      
//  cout << "at end " << in.tellg() << "\n";
}


void MfixData::ReadTimeValues(ifstream & in , int offset , int spxNum)
{
    in.clear();
    in.seekg( 3*512, ios::beg ); // first time

    float time;

    map<float,ULONGLONG> tmap2;

    while (in.read( (char*)&time,sizeof(float) ) )
    {
        SWAP_FLOAT(time);

        ULONGLONG pos = in.tellg();
//        cout << "------> " << time << "  " << pos << "\n";

        tmap2.insert( make_pair(time,pos - sizeof(float)) );
        times.insert(time);

        in.seekg(offset,ios::cur);
    }

    timeMap2[spxNum] = tmap2;
}


void MfixData::GetTimes()
{
    timeMap2.resize(nspx_use);

    for (int i=0; i<nspx_use; ++i)
    {
        int nvars = 0;

        string fname = run_name + ".SP" + ext[i];

        ifstream in(fname.c_str(),ios::binary);

        if (in) // file exists
        {

            switch (i+1)
            {

            case 1: nvars = 1; break;

            case 2: nvars = 2; break; 

            case 3: nvars = 3; break;

            case 4: nvars = 3*MMAX; break;

            case 5: nvars = MMAX; break;

            case 6:
                {
                    if (version_number <= 1.15)
                        nvars = 3;
                    else
                        nvars = MMAX + 1;

                    break;
                }

            case 7:
                {
                    nvars = NMAX[0];

                    for (int m=1; m<=MMAX; ++m) // bug
                    {
                        nvars += NMAX[m];
                    }

                    break;
                }


            case 8: nvars = MMAX; break;


            case 9: nvars = NScalar; break;

            case 10: nvars = nRR; break; 

            case 11:
                {
                    if (bKepsilon) nvars = 2;

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }
            }

            if (nvars > 0) 
            {
                int offset = 512-sizeof(float) + 512*(nvars*spx_records_per_timestep);
                ReadTimeValues( in , offset , i );
            }
        }
    }

    times.Sort(); // added SimpleSet

    nTimes = times.size();

}


void MfixData::CreateVariableNames()
{

    for (int i=0; i<nspx_use; ++i)
    {
        string fname = run_name + ".SP" + ext[i];

        ifstream in(fname.c_str(),ios::binary);

        if (in) // file exists
        {
            switch (i+1)
            {

            case 1:
                {
                    variable_names.push_back("EP_g");
                    varNameToSPX.insert(make_pair(string("EP_g"),1));
                    index_map.push_back(variable_names.size()-1);
                    break;
                }

            case 2:
                {
                    variable_names.push_back("P_g");
                    index_map.push_back(variable_names.size()-1);
                    variable_names.push_back("P_star");
                    index_map.push_back(variable_names.size()-1);
                    varNameToSPX.insert(make_pair(string("P_g"),2));
                    varNameToSPX.insert(make_pair(string("P_star"),2));
                    break;
                }

            case 3:
                {
                    variable_names.push_back("U_g");
                    index_map.push_back(variable_names.size()-1);
                    variable_names.push_back("V_g");
                    variable_names.push_back("W_g");
                    varNameToSPX.insert(make_pair(string("U_g"),3));
                    varNameToSPX.insert(make_pair(string("V_g"),3));
                    varNameToSPX.insert(make_pair(string("W_g"),3));
                    break;
                }


            case 4:
                {
                    string us = "U_s_";
                    string vs = "V_s_";
                    string ws = "W_s_";

                    for (int i=0; i<MMAX; ++i)
                    {
                        variable_names.push_back(us+ConvertIntToString(i+1));
                        index_map.push_back(variable_names.size()-1);
                        variable_names.push_back(vs+ConvertIntToString(i+1));
                        variable_names.push_back(ws+ConvertIntToString(i+1));
                        varNameToSPX.insert(make_pair(us+ConvertIntToString(i+1),4));
                        varNameToSPX.insert(make_pair(vs+ConvertIntToString(i+1),4));
                        varNameToSPX.insert(make_pair(ws+ConvertIntToString(i+1),4));

                    }

                    break;
                }

            case 5:
                {
                    string rops = "ROP_s_";

                    for (int i=0; i<MMAX; ++i)
                    {
                        variable_names.push_back(rops+ConvertIntToString(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(rops+ConvertIntToString(i+1) , 5));
                    }

                    break;
                }

            case 6:
                {
                    variable_names.push_back("T_g");
                    index_map.push_back(variable_names.size()-1);
                    varNameToSPX.insert(make_pair(string("T_g"),6));

                    if (version_number <= 1.15)
                    {
                        variable_names.push_back("T_s_1");
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(string("T_s_1"), 6));

                        if (MMAX > 1) 
                            variable_names.push_back("T_s_2");
                        else
                            variable_names.push_back("T_s_2_not_used");

                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(string("T_s_2"),6));
                    }
                    else
                    {
                        for (int i=0; i<MMAX; ++i)
                        {
                            string ts = "T_s_";

                            for (int i=0; i<MMAX; ++i)
                            {
                                variable_names.push_back(ts+ConvertIntToString(i+1));
                                index_map.push_back(variable_names.size()-1);
                                varNameToSPX.insert(make_pair(ts+ConvertIntToString(i+1) , 6));
                            }
                        }
                    }

                    break;
                }

            case 7:
                {
                    string var = "X_g_";

                    for (int i=0; i<NMAX[0]; ++i)
                    {
                        variable_names.push_back(var+ConvertIntToString(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+ConvertIntToString(i+1) , 7));
                    }

                    var = "X_s_";

                    for (int m=0; m<MMAX; ++m)
                    {
                        for (int i=0; i<NMAX[m]; ++i)
                        {
                            string name = var + ConvertIntToString(m+1) + "_"
                                              + ConvertIntToString(i+1);

                            variable_names.push_back(name);
                            index_map.push_back(variable_names.size()-1);
                            varNameToSPX.insert(make_pair(name, 7));
                        }
                    }


                    break;
                }


            case 8:
                {
                    string var = "Theta_m_";

                    for (int i=0; i<MMAX; ++i)
                    {
                        variable_names.push_back(var+ConvertIntToString(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+ConvertIntToString(i+1),8));
                    }

                    break;
                }


            case 9:
                {
                    string var = "Scalar_";

                    for (int i=0; i<NScalar; ++i)
                    {
                        variable_names.push_back(var+ConvertIntToString(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+ConvertIntToString(i+1), 9));
                    }

                    break;
                }


            case 10:
                {
                    string var = "RRates_";

                    for (int i=0; i<nRR; ++i)
                    {
                        variable_names.push_back(var+ConvertIntToString(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+ConvertIntToString(i+1), 10));
                    }

                    break;
                }

            case 11:
                {
                    if (bKepsilon)
                    {
                        variable_names.push_back("k_turb_g");
                        index_map.push_back(variable_names.size()-1);
                        variable_names.push_back("e_turb_g");
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(string("k_turb_g"),11));
                        varNameToSPX.insert(make_pair(string("e_turb_g"),11));
                    }

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }


            }


        }
    }

}


void MfixData::GetVariableAtTimestep(int vari , int tstep)
{
    // This routine opens and closes the file for each request.
    // Maybe keep all SPX files open, and just perform relative
    // moves to get to the correct location in the file

    // get filename that vaiable # vari is located in

    // assumptions : there are <10 solid phases,
    // <10 scalars and <10 ReactionRates (need to change this)

    string & vname = variable_names[vari];

    int spx = varNameToSPX[vname];

    string fname = run_name + ".SP" + ext[spx-1];

    // spx_timeMap2 ... map of times and file position for the
    //                 specified SPX file

    map<float,ULONGLONG> & spx_timeMap2 = timeMap2[spx-1];

    // find where the requested time  , times[tstep] , is located in the map

    map<float,ULONGLONG>::iterator mit = spx_timeMap2.lower_bound(times[tstep]);

    if (mit == spx_timeMap2.end()) // time is greater than last time in file
    {
        mit = --spx_timeMap2.end(); // point to last element
    }
    else
    {
        // mit either points to correct time if the time exists
        // in the file , or the next larger time. For now, I want
        // to return the first time before the requested time
        // if it is not in the file 

        // could do time interpolation instead
        // could use the closest time instead

        if (mit->first > times[tstep] && mit->second != 0) --mit;
    }

    int file_position = mit->second;

    int nskip = 0;  // number of variables to skip to get to desired variable
                    // example : request : U_g ... nskip = 0
                    // example : request : V_g ... nskip = 1
                    // example : request : W_g ... nskip = 2

    switch (spx)
    {

    case 1: nskip = 0; break;

    case 2:
        {
            if (vname == "P_g")
                nskip = 0;
            else
                nskip = 1;

            break;
        }

    case 3:
        {
            if (vname == "U_g")
                nskip = 0;
            else if (vname == "V_g")
                nskip = 1;
            else
                nskip = 2;

            break;
        }

    case 4:
        {
            // form of variable : "A_s_b"
            // where "A" = U / V / W
            // and "b" = solid phase index ( 1 - based )

            // examples : U_s_1 , W_s_2

            int m = ConvertCharToInt(vname[ vname.size()-1 ]);

            nskip = 3*(m-1);
            
            if (vname[0] == 'V') ++nskip;
            if (vname[0] == 'W') nskip += 2;

            break;
        }

    case 5:
        {
            // form of variable : "ROP_s_b"
            // "b" = solid phase index ( 1 - based )

            int m = ConvertCharToInt(vname[ vname.size()-1 ]);
            nskip = m - 1;

            break;
        }

    case 6:
        {
            // form of variable : "T_g" or "T_s_b"
            // "b" = solid phase index ( 1 - based )

            if (vname == "T_g")
            {
                nskip = 0;
            }
            else
            {
                int m = ConvertCharToInt(vname[ vname.size()-1 ]);
                nskip = m; // 1 + (m-1)
            }

            break;
        }

    case 7:
        {
            // form of variable : "X_g_a" or "X_s_b_c"
            // "a" = gas species index ( 1 - based )
            // "b" = solid phase index ( 1 - based )
            // "c" = solid species index ( 1 -based )

            if (vname.substr(0,3) == "X_g")
            {
                int ng = ConvertCharToInt(vname[ vname.size()-1 ]);
                nskip = ng - 1;
            }
            else
            {
                int ns = ConvertCharToInt(vname[ vname.size()-1 ]);

                string::size_type pos1 = vname.find_first_of('_');
                string::size_type pos2 = vname.find_last_of('_');

                string str_m = vname.substr( pos1+3 , pos2-pos1-3 );
                int m = ConvertStringToInt(str_m);

                nskip = NMAX[0];

                for (int mm=1; mm<m-1; ++mm)
                {
                    nskip += NMAX[mm];
                }

                nskip += ns - 1;
            }

            break;
        }

    case 8:
        {
            // form of variable : "Theta_m_b"
            // "b" = solid phase index ( 1 - based )

            int m = ConvertCharToInt(vname[ vname.size()-1 ]);
            nskip = m - 1;

            break;
        }

    case 9:
        {
            // form of variable : "Scalar_s"
            // "s" = scalar index ( 1 - based )

            int m = ConvertCharToInt(vname[ vname.size()-1 ]);
            nskip = m - 1;

            break;
        }

    case 10:
        {
            // form of variable : "RRates_r"
            // "r" = reaction rate index ( 1 - based )

            int m = ConvertCharToInt(vname[ vname.size()-1 ]);
            nskip = m - 1;

            break;
        }

    case 11:
        {
            if (vname == "k_turb_g")
                nskip = 0;
            else
                nskip = 1;

            break;
        }
    } // switch

    // get the requested variable/time data

    ULONGLONG nBytesSkip = mit->second;

    nBytesSkip += 512; // record with time in it

    nBytesSkip += nskip*spx_records_per_timestep*512;

    ifstream in(fname.c_str(),ios::binary);

    in.seekg(nBytesSkip,ios::beg);

    IN_BIN_512R (in,&var[0],ijkmax2);
}




void MfixData::Split_SP1(int SPX_file)
{
//    for (int i=0; i<nspx_use; ++i)
    for (int i=SPX_file; i<SPX_file+1; ++i)   // just do one file
    {
        int nvars = 0;

        string fname = run_name + ".SP" + ext[i];

	string fname_a = run_name + ".SP" + ext[i] + "a";
	string fname_b = run_name + ".SP" + ext[i] + "b";

        ifstream in(fname.c_str(),ios::binary);

        if (in) // file exists
        {

            switch (i+1)
            {

            case 1: nvars = 1; break;

            case 2: nvars = 2; break; 

            case 3: nvars = 3; break;

            case 4: nvars = 3*MMAX; break;
   int NVARS_SPX_file(int SPX_file);

            case 5: nvars = MMAX; break;

            case 6:
                {
                    if (version_number <= 1.15)
                        nvars = 3;
                    else
                        nvars = MMAX + 1;

                    break;
                }

            case 7:
                {
                    nvars = NMAX[0];

                    for (int m=1; m<=MMAX; ++m) // bug
                    {
                        nvars += NMAX[m];
                    }

                    break;
                }


            case 8: nvars = MMAX; break;


            case 9: nvars = NScalar; break;

            case 10: nvars = nRR; break; 

            case 11:
                {
                    if (bKepsilon) nvars = 2;

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }
            }

            if (nvars > 0) 
            {
                int offset = 512-sizeof(float) + 512*(nvars*spx_records_per_timestep);
            //    ReadTimeValues( in , offset , i );
                ReadWriteValues( in , i , nvars , fname_a , fname_b);
            }
        }
    }
}


void MfixData::ReadWriteValues(ifstream & in , int spxNum , int nvars ,
             const string & fname_a , const string & fname_b)
{
    ofstream out1(fname_a.c_str(),ios::binary);
    ofstream out2(fname_b.c_str(),ios::binary);

    unsigned int c = 3;

    in.clear();
    in.seekg( 0, ios::beg ); 
   int NVARS_SPX_file(int SPX_file);

    char buf[512];

    // read & write the first 3 records (header)
    for (int i=0; i<3; ++i)
    {
       in.read( buf , 512 );
       out1.write( buf , 512 );
       out2.write( buf , 512 );
    }

    int rec1 = 4;
    int rec2 = 4;

    bool bSwitched = false;

    float time , time_split , time_start;

    cout << "enter the time that should be used to split the file > ";
    cin  >> time_split;

    cout << "enter the time to start writing to the file > ";
    cin  >> time_start;

    while (in.read( (char*)&time,sizeof(float) ) )
    {

        SWAP_FLOAT(time);

        float time_save = time;

        if (time <= time_split)
        {
           if (time >= time_start) ++rec1;
        }
        else
        {
           if (time >= time_start) ++rec2;
        }

        ++c;

        cout << "-processing time = " << time << "   (record # " << c << ")\n";


        in.read( buf , 512-sizeof(float) );
        SWAP_FLOAT(time);

        if (time_save <= time_split)
        {
           if (time_save >= time_start)
           {
              out1.write( (char*)(&time),sizeof(float) );
              out1.write( buf , 512-sizeof(float) );
           }
        }
        else
        {
           if (!bSwitched)
           {
              cout << "switching to second file\n";
              bSwitched = true;
              out1.close();

              cout << fname_a << " ... updating record 3 ... n = " << rec1 << "\n";
              SWAP_INT(rec1);
              fstream out(fname_a.c_str(),ios::binary|ios::in|ios::out);
              out.seekp(512*2,ios::beg);
              out.write( (char*)(&rec1) , 4);
              out.close();
              SWAP_INT(rec1);

           }

           if (time_save >= time_start)
           {
              out2.write( (char*)(&time),sizeof(float) );
              out2.write( buf , 512-sizeof(float) );
           }
        }
 
        if (time_save >= time_start)
        {
		cout << " in write mode\n";
        	for (int i=0; i<nvars*spx_records_per_timestep; ++i)
        	{
           		++c;
           		in.read( buf , 512 );
           		if (time_save <= time_split)

           		{
             			++rec1;
             			out1.write( buf , 512 );
           		}
           		else
           		{
             			++rec2;
             			out2.write( buf , 512 );
           		}
        	}
        }
        else
        {
             in.seekg( 512*nvars*spx_records_per_timestep , ios::cur);
             c += nvars*spx_records_per_timestep;
        }
        cout << "rec1/2 = " << rec1 << " : " << rec2 << "\n";
    }

    if (!bSwitched)
    {
        out1.close();

        cout << fname_a << " ... updating record 3 ... n = " << rec1 << "\n";
        SWAP_INT(rec1);
        fstream out(fname_a.c_str(),ios::binary|ios::in|ios::out);
        out.seekp(512*2,ios::beg);
        out.write( (char*)(&rec1) , 4);
        out.close();
        SWAP_INT(rec1);
    }

    out2.close();

    cout << fname_b << " ... updating record 3 ... n = " << rec2 << "\n";
    fstream out(fname_b.c_str(),ios::binary|ios::in|ios::out);
    out.seekp(512*2,ios::beg);
    SWAP_INT(rec2);
    out.write( (char*)(&rec2) , 4);
    out.close();

   
   ifstream in1(fname_a.c_str(),ios::binary);
   in1.seekg(2*512,ios::beg);
   in1.read( (char*)(&rec1),4);
   SWAP_INT(rec1);
   cout << "first file ... check ... rec4 = " << rec1 << "\n";

   ifstream in2(fname_b.c_str(),ios::binary);
   in2.seekg(2*512,ios::beg);
   in2.read( (char*)(&rec2),4);
   SWAP_INT(rec2);
   cout << "second file ... check ... rec4 = " << rec2 << "\n";
}

int MfixData::NVARS_SPX_file(int SPX_file)
{
   int nvars = -1;

   switch (SPX_file+1)
  {

            case 1: nvars = 1; break;

            case 2: nvars = 2; break; 

            case 3: nvars = 3; break;

            case 4: nvars = 3*MMAX; break;

            case 5: nvars = MMAX; break;

            case 6:
                {
                    if (version_number <= 1.15)
                        nvars = 3;
                    else
                        nvars = MMAX + 1;

                    break;
                }

            case 7:
                {
                    nvars = NMAX[0];

                    for (int m=1; m<=MMAX; ++m) // bug
                    {
                        nvars += NMAX[m];
                    }

                    break;
                }


            case 8: nvars = MMAX; break;


            case 9: nvars = NScalar; break;

            case 10: nvars = nRR; break; 

            case 11:
                {
                    if (bKepsilon) nvars = 2;

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << SPX_file+1 << "\n";
                    break;
                }
   }
   return nvars;
}

float MfixData::GetLastTime(int SPX_file)
{
	string fname = run_name + ".SP" + ext[SPX_file];
        ifstream in(fname.c_str(),ios::binary);

        long pos = 512 + 512*NVARS_SPX_file(SPX_file)*spx_records_per_timestep;

        in.seekg(0,ios::end);

        in.seekg(-pos,ios::cur);

        float t;
        GetFloat(in,t);

	return t;
}

float MfixData::GetLastTime(fstream & in , int nvars)
{
        long pos = 512 + 512*nvars*spx_records_per_timestep;

	in.clear();

        in.seekg(0,ios::end);

        in.seekg(-pos,ios::cur);

        float t;
        GetFloat(in,t);

	return t;
}

int MfixData::GetRec3_value(fstream & in)
{
	in.clear();

	in.seekg(512*2,ios::beg);

        int rec3;

        GetInt(in,rec3);

	return rec3;
}



void MfixData::CombineSPX(int SPX_file , MfixData & data2)
{
	if ( NVARS_SPX_file(SPX_file)  != data2.NVARS_SPX_file(SPX_file) )
	{
		cout << "number of variables in the SPx file are not equal\n\n";
		return;
	}

	string fname = run_name + ".SP" + ext[SPX_file];
        fstream in(fname.c_str(),ios::binary | ios::in | ios::out);

        if (!in)
	{
		cout << "could not open file : " << fname << "\n";
		return;
	}

	int nvars =  NVARS_SPX_file(SPX_file);

	int   rec3      = GetRec3_value(in);
	float last_time = GetLastTime(in,nvars);

        cout << "\n";
  //      cout << "rec3      = " << rec3      << "\n";
        cout << "last time (first file) = " << last_time << "\n";

	// open second file and get to first time record
	string fname2 = data2.run_name + ".SP" + ext[SPX_file];
        ifstream in2(fname2.c_str(),ios::binary | ios::in | ios::out);

	cout << fname  << "\n";
	cout << fname2 << "\n\n\n";


	// skip first time in second file ... it is a repeat
	// of the last time in the first file

        ULONGLONG pos  = 3*512;  				// header
                  pos += 512;    				// first time
                  pos += 512*nvars*spx_records_per_timestep;	// data for first time

        in2.seekg(pos,ios::beg);

	// process the second file ... appending to the first

        in.seekp(0,ios::end);

	char  buf[512];
	float time2 , time_file;

        while (in2.read(buf,512))
        {
		time2 = *(reinterpret_cast<float*>(buf));
                SWAP_FLOAT(time2);
                time_file = time2 + last_time;
		cout << "processing time = " << setw(15) << time_file 
                                             << " = " 
                                             << setw(15) << last_time
                                             << " + "
                                             << time2
                                             << "\n";

 		SWAP_FLOAT(time_file);

		++rec3;
		in.write(reinterpret_cast<char*>(&time_file),sizeof(float));
		in.write(buf+4,512-sizeof(float));

		for (int i=0; i<nvars*spx_records_per_timestep; ++i)
		{
			in2.read(buf,512);
			in.write(buf,512);
                        ++rec3;
		}

        }

	// write the "next record" value in record 3

        cout << "updating record 3 value ...\n";

	in.seekp(2*512,ios::beg);
	SWAP_INT(rec3);
	in.write( reinterpret_cast<char*>(&rec3),sizeof(int) );
}


