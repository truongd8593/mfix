#include <QtGui>

#include <string>

#include "ani_mfix.h"
#include "MfixData.h"
#include "data.h"

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
	 string run_name = it->toStdString();
	 run_name = run_name.substr( 0 , run_name.size() - 4 ); // remove the .RES
	 
         Data::mfix.SetName(run_name.c_str());
	 Data::mfix.ReadRes0();
         Data::mfix.CreateVariableNames();
         Data::mfix.GetTimes();
	 
	 Data::SetVariableRanges();
	 Data::vPlotWindowInfo.resize(9);
	 
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

  
