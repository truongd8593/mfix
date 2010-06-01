#include <QtGui>

#include <string>

#include "ani_mfix.h"
#include "MfixData.h"
#include "data.h"

extern void Read_cut_plane_info(const char * fname);
extern void set_geom(double * ddx , const int & nx ,
					double * ddy , const int & ny ,
					double * ddz , const int & nz );


int main(int argc , char * argv[])
{
   using namespace std;
   
   Data::init();
   
   QApplication app(argc,argv);

   QFileDialog fd(0, QString("chose MFIX restart file") , QString::null , QString("MFIX (*.RES)") );
   fd.show();
   
   if (fd.exec() )
   {
      QStringList filenames = fd.selectedFiles();
      
      QStringList::iterator it = filenames.begin();
      
      if (it != filenames.end())
      {
	 QDir qdir      = fd.directory();
         QString qs_dir = qdir.absolutePath();
         Data::sCWD     = qs_dir.toStdString();
         const size_t dir_size = Data::sCWD.size();
         if ( dir_size > 2)
         {
            if (Data::sCWD[ dir_size - 1] != '/') Data::sCWD.push_back('/');
         }
 
	 string run_name = it->toStdString();
	 run_name = run_name.substr( 0 , run_name.size() - 4 ); // remove the .RES
	 
         Data::mfix.SetName(run_name.c_str());
	 Data::mfix.ReadRes0();
         Data::mfix.CreateVariableNames();
         Data::mfix.GetTimes();
	 
	 Data::SetVariableRanges();
	 Data::vPlotWindowInfo.resize(9);
	 
	 
	 set_geom (
	            &Data::mfix.DX[0] , Data::mfix.imax2 ,
	            &Data::mfix.DY[0] , Data::mfix.jmax2 ,
	            &Data::mfix.DZ[0] , Data::mfix.kmax2
                  );

         string vtk_name = Data::sCWD + "ani_cutcell.vtk";
	 Read_cut_plane_info(vtk_name.c_str());
	 
         AniMfix * client = new AniMfix;
         client->setAttribute(Qt::WA_DeleteOnClose);
         client->show();
         int ret_val = app.exec();
       //  delete client;
         return ret_val;
	 
      }
   }
   
   return 0;
}

  
