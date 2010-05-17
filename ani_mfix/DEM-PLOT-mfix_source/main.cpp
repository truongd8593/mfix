#include <QtGui>

#include <string>

#include "ani_mfix.h"
#include "data.h"
#include "MfixData.h"

int main(int argc , char * argv[])
{
   using namespace std;
      
   Data::init();
   
   if (argc > 1) Data::bFake_diameter = true;
   
   QApplication app(argc,argv);

   QFileDialog fd(0, QString("chose MFIX restart file") , QString::null , QString("MFIX (*.RES)") );
   fd.show();
   
   if (fd.exec() )
   {
      QStringList filenames = fd.selectedFiles();
      
      QStringList::iterator it = filenames.begin();
      
      if (it != filenames.end())
      {
	 Data::run_name = it->toStdString();
	 Data::run_name = Data::run_name.substr( 0 , Data::run_name.size() - 4 ); // remove the .RES
	 
	 string mfix_name = Data::run_name;
	 Data::run_name += "_DES";
	 
	 Data::mfix.SetName(mfix_name.c_str());
	 Data::mfix.ReadRes0();
	 Data::mfix_xlength = Data::mfix.xlength;
	 Data::mfix_ylength = Data::mfix.ylength;
	 Data::mfix_zlength = Data::mfix.zlength;
         Data::mfix_imax2   = Data::mfix.imax2;
         Data::mfix_jmax2   = Data::mfix.jmax2;
         Data::mfix_kmax2   = Data::mfix.kmax2;

	 if (Data::mfix.coordinates[1] == 'Y' || Data::mfix.coordinates[1]=='y' )	
	    Data::mfix_bCartesian = false;
	 else
	    Data::mfix_bCartesian = true;

	 if ( Data::mfix.kmax2>1 && (Data::mfix.coordinates[1] == 'Y' || Data::mfix.coordinates[1]=='y') )
	 {
	   Data::xmin = -Data::mfix.xlength;
	 }
	 else
	 {
	   Data::xmin = 0;
	 }
 
	 Data::xmax = Data::mfix.xlength;
	 Data::ymin = 0;
	 Data::ymax = Data::mfix.ylength;

	 Data::GetTimes();
	 Data::cur_ts = 0;
	 Data::ReadNext_XML_file();

	 	 
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

  
