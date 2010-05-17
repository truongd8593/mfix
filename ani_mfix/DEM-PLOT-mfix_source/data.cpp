#include <QtOpenGL>
#include <QtXml/QXmlStreamReader>
#include <QXmlStreamAttributes>

#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>
#include <set>
#include <iomanip>
#include <cctype>

#include "data.h"

bool Data::bPaused        = true;
bool Data::bBusy          = false;
bool Data::b3D_bubbles    = true;
bool Data::bAdd_to_vpos   = true;
bool Data::bTrace         = false;
bool Data::bFake_diameter = false;
bool Data::bImage         = false;


MfixData Data::mfix; //

int  Data::gl_width = 10;
int  Data::gl_height = 10;
int  Data::npixels_x = 10;
int  Data::npixels_y = 10;
int  Data::cur_ts = 0;
int  Data::color_method = 0;

int  Data::sphere_resolution = 7;

std::string Data::run_name;

std::vector<Particle> Data::vP;
std::vector<Particle> Data::vPos;
std::vector<double>   Data::vTimes;
std::vector<int>      Data::vColors_position;
std::vector<int>      Data::vColors_diameter;

double Data::vel_max = 0.00001;

float Data::xmin = 0;
float Data::xmax = 1;

float Data::ymin = 0;
float Data::ymax = 1;

float Data::scale = 1.0f;
float Data::xtrans = 0;
float Data::ytrans = 0;
int   Data::xrot   = 0;
int   Data::yrot   = 0;

double Data::bAdd_time = -1;

double Data::time = 0;

MyFont Data::myfont;

float Data::mfix_xlength    = 1;
float Data::mfix_ylength    = 1;
float Data::mfix_zlength    = 1;
bool  Data::mfix_bCartesian = 1;
int   Data::mfix_imax2      = 1;
int   Data::mfix_jmax2      = 1;
int   Data::mfix_kmax2      = 1;



namespace
{
  const int title_height = 75;
  const int time_height  = 75;
  
  /*
  float rra[] = { 1.00 , 0.95 , 0.90 , 0.90 , 0.90 , 0.90 , 
       0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.75 , 0.65 , 0.55 , 0.50 , 0.40 , 0.40 , 0.40 , 
       0.35 , 0.35 , 0.30 , 0.30 , 0.30 , 0.30 , 0.20 , 
       0.20 , 0.20 , 0.10 , 0.00 , 0.00 , 1.00 }; 
 
   float gga[] = { 0.00 , 0.10 , 0.15 , 0.40 , 0.55 , 0.65 , 
       0.75 , 0.80 , 0.85 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.85 , 0.80 , 0.75 , 0.65 , 
       0.55 , 0.40 , 0.15 , 0.10 , 0.00 , 1.0 };
 
   float bba[] = { 0.00 , 0.00 , 0.10 , 0.20 , 0.20 , 0.20 , 
       0.30 , 0.30 , 0.30 , 0.30 , 0.35 , 0.35 , 0.40 , 
       0.40 , 0.40 , 0.50 , 0.55 , 0.65 , 0.75 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,  
       0.90 , 0.90 , 0.90 , 0.95 , 1.00 , 1.0};
       */

    std::set<float> sDiameters;
}


void Data::init()
{
  cur_ts = 0;
  b3D_bubbles = true;
  bAdd_to_vpos = true;
  bTrace       = false;
  scale = 1.0f;
  xtrans = 0;
  ytrans = 0;
  xrot   = 0;
  yrot   = 0;
  time = 0;
  bAdd_time = -1;
  color_method = 0;
  vel_max = 0.00001;
  bFake_diameter = false;
  bImage = false;
  myfont.Process();

  mfix_xlength = 1;
  mfix_ylength = 1;
  mfix_zlength = 1;
  mfix_bCartesian = true;
  mfix_imax2   = 1;
  mfix_jmax2   = 1;
  mfix_kmax2   = 1;
  
  sphere_resolution = 7;
}

struct Lower
{
  char operator() (char c) { return ::tolower(c); }
};

void Data::ReadNext_XML_file()
{
    using namespace std;
    
    vP.clear();
  
    std::stringstream ss;
    ss << Data::run_name << "_" << setw(5) << setfill('0') << Data::cur_ts << ".vtp";
    string fname = ss.str();
   
    QFile* file = new QFile(fname.c_str());    
    
    if (!file->open(QIODevice::ReadOnly | QIODevice::Text)) 
    {
      cur_ts = 0;
      vPos.clear();
      delete file;
      return;
    }    
    
    QXmlStreamReader xml(file);
    
    while(!xml.atEnd() && !xml.hasError()) 
    {
	QXmlStreamReader::TokenType token = xml.readNext();  // read next element
	
	// If token is just StartDocument, we'll go to next
	    
	if(token == QXmlStreamReader::StartDocument) 
	{
	  continue;
	}
      
      // processing instructuion (Time= in our case)
      if (token == QXmlStreamReader::ProcessingInstruction) 
      {
        std::string starget = xml.processingInstructionTarget().toString().toStdString();
	transform(starget.begin(),starget.end(),starget.begin(),Lower());
	
	if (starget == "time")
	{
	    std::string s = xml.processingInstructionData().toString().toStdString();
	    std::replace(s.begin(),s.end(),'=',' ');
	    Data::time = Convert<std::string,double>(s);
	}
	

      }
      
      // If token is StartElement, we'll see if we can read it.
      if(token == QXmlStreamReader::StartElement) 
      {
	
	if (xml.name() == "Piece")
	{
	  QXmlStreamAttributes qa = xml.attributes();
	  
	  for (int i=0; i<qa.count(); ++i)
	  {
	      if (qa[i].name().toString() == "NumberOfPoints")
	      {
		int n = Convert<string,int>( qa[i].value().toString().toStdString() );
		vP.resize(n);
		vColors_position.resize(n,1);
		vColors_diameter.resize(n,1);
	      }
	  }
	  
	}
	// If it's named persons, we'll go to the next.
	QString s = xml.name().toString();
	
	if(xml.name() == "DataArray") 
	{
	    QXmlStreamAttributes qa = xml.attributes();
 
	    bool bDiameter  = false;
	    bool bVelocity  = false;
	    bool bPosition  = false;
	    
	    for (int i=0; i<qa.count(); ++i)
	    {
	      if (qa[i].name() == "Name" || qa[i].name()=="NAME"  )
	      {
		  if (qa[i].value() == "Diameter") bDiameter = true;
		  if (qa[i].value() == "Velocity") bVelocity = true;
		  if (qa[i].value() == "Position") bPosition = true;
	      }
	    }
	    	    
	    QString s = xml.readElementText();
	
	    std::stringstream ss( s.toStdString() );
	    std::string line;
	   
	    size_t c = 0;
	    while (getline(ss,line) &&  c<vP.size())
	    {
		double v1,v2,v3;
	      
		if (line.size() > 6)
		{
		    stringstream sss(line);
		    
		    sss >> v1 >> v2 >> v3;
		    
		    if (bVelocity)
		    {
			vP[c].u = v1;
			vP[c].v = v2;
			vP[c].w = v3;
                      //  if (cur_ts == 0)
			{
			  double vel = v1*v1 + v2*v2 + v3*v3;
			  if (vel > vel_max) vel_max = vel;
			}
		    }
		    else if (bPosition)
		    {
			vP[c].x = v1;
			vP[c].y = v2;
			vP[c].z = v3;

                        if (cur_ts == 0)
			{
			    if (v1 <= (xmin+xmax)/2)
				vColors_position[c] = 0;
			    else
				vColors_position[c] = 1;
			}
		    }
		    else if (bDiameter)
		    {
		      
		      if (bFake_diameter)
		      {
			  if (c<vP.size()/4) 
			    v1 += 0.01; // just for testing
			  else if (c<vP.size()/2)
			    v1 += 0.02; // just for testing
			  else if (c<3*vP.size()/4)
			    v1 += 0.03; // just for testing
			  else 
			    v1 += 0.04; // just for testing
		      }

		      vP[c].d = v1;
		      if (cur_ts == 0)
		      {
			  sDiameters.insert(v1);
		      }
		    }
		    
		    ++c;
		}
	    }
	    continue;
	}
//		    /* If it's named person, we'll dig the information from there.*/
	if(xml.name() == "DataArray") 
	{
	}
      }
      else
      {
      }
    }
    /* Error handling. */
    if(xml.hasError()) 
    {
      std::cout << "an error occurred\n";
    }
    xml.clear();
    
    delete file;

    if (cur_ts == 0)
    {
	for (size_t i=0; i<vColors_diameter.size(); ++i)
	{
	    set<float>::iterator it = sDiameters.find(vP[i].d);
	    size_t dis = std::distance(sDiameters.begin(),it);
	    vColors_diameter[i] = dis;
	}
    }
}

void Data::GetTimes()
{
    using namespace std;
    
    
    bool bGo = true;
    
    cur_ts = 0;
    
    while (bGo)
    {
	  std::stringstream ss;
	  ss << Data::run_name << "_" << setw(5) << setfill('0') << Data::cur_ts << ".vtp";
	  string fname = ss.str();
	
	  QFile* file = new QFile(fname.c_str());    
	  
	  if (!file->open(QIODevice::ReadOnly | QIODevice::Text)) 
	  {
	    cur_ts = 0;
	    delete file;
	    return;
	  }    
	  
	  QXmlStreamReader xml(file);
	  
	  while(!xml.atEnd() && !xml.hasError()) 
	  {
	      QXmlStreamReader::TokenType token = xml.readNext();  // read next element
	      
	      // If token is just StartDocument, we'll go to next
		  
	    
	    // processing instructuion (Time= in our case)
	    if (token == QXmlStreamReader::ProcessingInstruction) 
	    {
	      std::string starget = xml.processingInstructionTarget().toString().toStdString();
	      transform(starget.begin(),starget.end(),starget.begin(),Lower());
	      
	      if (starget == "time")
	      {
		  std::string s = xml.processingInstructionData().toString().toStdString();
		  std::replace(s.begin(),s.end(),'=',' ');
		  Data::time = Convert<std::string,double>(s);
		  vTimes.push_back( Data::time );
		  ++Data::cur_ts;
		  break;
	      }
	    }
	    
	  }
	  /* Error handling. */
	  if(xml.hasError()) 
	  {
	    std::cout << "an error occurred\n";
	  }
	  xml.clear();
	  
	  delete file;
    }
}

