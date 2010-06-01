#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <map>
#include <set>
#include <cmath>

#include <QtOpenGL>

#include "data.h"


using namespace std;

//#define for if (0); else for

//ofstream pan("zzzzzz.txt");


struct CutPoint
{
	double x,y,z;

	CutPoint( double xx=0 , double yy=0 , double zz=0) : x(xx) , y(yy) , z(zz) {}
};


struct Vector
{
	double x,y,z;

	Vector() {}

	Vector(const CutPoint & rhs) : x(rhs.x) , y(rhs.y) , z(rhs.z) {}

};

CutPoint operator + (const CutPoint & p1 , const CutPoint & p2)
{
	CutPoint p;
	p.x = p1.x + p2.x;
	p.y = p1.y + p2.y;
	p.z = p1.z + p2.z;

	return p;
}

CutPoint operator - (const CutPoint & p1 , const CutPoint & p2)
{
	CutPoint p;
	p.x = p1.x - p2.x;
	p.y = p1.y - p2.y;
	p.z = p1.z - p2.z;

	return p;
}

Vector operator * (double sc, const Vector & rhs)
{
	Vector v;
	v.x = sc * rhs.x;
	v.y = sc * rhs.y;
	v.z = sc * rhs.z;

	return v;
}

Vector operator + (const Vector & p1 , const Vector & p2)
{
	Vector p;
	p.x = p1.x + p2.x;
	p.y = p1.y + p2.y;
	p.z = p1.z + p2.z;

	return p;
}

Vector operator - (const Vector & p1 , const Vector & p2)
{
	Vector p;
	p.x = p1.x - p2.x;
	p.y = p1.y - p2.y;
	p.z = p1.z - p2.z;

	return p;
}


template <class charT, class traits>
inline std::basic_istream<charT,traits>&
ignoreLine (std::basic_istream<charT,traits>& strm)
{
    strm.ignore(100000,strm.widen('\n'));

    return strm;
}






extern "C"
{


	void plot_ij_slice_(const int & kuse , float * ra , float * ga , float * ba ,
		                int * color_array , const int & fp , 
                                const int * flag , int * status);

}


struct CutPolygon
{
	vector<int> vIndex;
};


struct CutMfix
{
	vector<CutPoint> vCenterOfMass;

	void CalculateCenterOfMass(const vector<CutPoint>   & vPoints ,
		                       const vector<CutPolygon> & vPolygons);
};

struct IJ
{
	int i,j;

	IJ(int ii=-1,int jj=-1) : i(ii) , j(jj) {}

	bool operator < (const IJ & rhs) const
	{
		if (i<rhs.i) return true;
		if (i>rhs.i) return false;
		return j<rhs.j;
	}
};


istream & operator >> (istream & in , CutPoint & cp)
{
	in >> cp.x >> cp.y >> cp.z >> ignoreLine;

	return in;
}

ostream & operator << (ostream & out , const CutPoint & cp)
{
	out << "(" << cp.x << " , " <<  cp.y << " , " <<  cp.z << ")";

	return out;
}

istream & operator >> (istream & in , CutPolygon & cp)
{
	int n;

	in >> n;

	cp.vIndex.resize(n);

	for (int i=0; i<n; ++i)
	{
		in >> cp.vIndex[i];
	}

	in >> ignoreLine;

	return in;
}


namespace
{
	vector<CutPoint>   vPoints;
	vector<CutPolygon> vPolygons;

	CutMfix            cut_mfix;

	vector<double> dx , dy , dz;

	map<IJ,int>        map_IJ_to_vpolygons;

	set<int>           set_k;

	void SwapEndian(double & x)
	{
    	char * p = reinterpret_cast<char*>(&x);

    		swap(p[0],p[7]);
    		swap(p[1],p[6]);
    		swap(p[2],p[5]);
    		swap(p[3],p[4]);
	}

	void SwapEndian(int & x)
	{
    		char * p = reinterpret_cast<char*>(&x);

    		swap(p[0],p[3]);
    		swap(p[1],p[2]);
	}

        bool bNo_vtk = true;


}

void Read_cut_plane_info(const char * fname)
{

    ifstream in(fname,ios::binary);

    string geom_name = Data::sCWD + "cutcells_geom.txt";
    ofstream out(geom_name.c_str());

    if (!in) 
    {
        cout << fname << " : cutfile geometry not found ... using normal processing\n";
        return;
    }

    string line;

    getline(in,line);
    getline(in,line);
    getline(in,line);
    getline(in,line);
    getline(in,line);
 
    double x,y,z;

    int c = 0;

    while (in.peek() != '\n')
    {
        ++c;
        in.read(reinterpret_cast<char*>(&x),sizeof(double));
        in.read(reinterpret_cast<char*>(&y),sizeof(double));
        in.read(reinterpret_cast<char*>(&z),sizeof(double));

        SwapEndian(x);
        SwapEndian(y);
        SwapEndian(z);

	CutPoint cp;
	cp.x = x;
	cp.y = y;
	cp.z = z;
	vPoints.push_back(cp);


    }

    out << vPoints.size() << " = numbers of points that follow\n";
    for (size_t i=0; i<vPoints.size(); ++i)
    {
       out << vPoints[i].x << " " << vPoints[i].y << " " << vPoints[i].z << "\n";
    }

 //   cout << "vPoints.size() = " << vPoints.size() << "\n";

    getline(in,line);  // get rid of the '\n'

    getline(in,line);
//    cout << line << "\n";

    int n,index;

    c = 0;
    
    while (in && in.peek() != '\n')
    {
        ++c;
        in.read(reinterpret_cast<char*>(&n),sizeof(int));

        SwapEndian(n);

        CutPolygon cp;
	cp.vIndex.reserve(n);


        for (int i=0; i<n; ++i)
        {
            in.read(reinterpret_cast<char*>(&index),sizeof(int));

            SwapEndian(index);
            cp.vIndex.push_back(index);
        }

	vPolygons.push_back(cp);
    }

    out << vPolygons.size() << " = numbers of polygons that follow\n";
    for (size_t i=0; i<vPolygons.size(); ++i)
    {
       out << vPolygons[i].vIndex.size() << " = number of vertices for polygon # " << i << "\n";
       for (size_t j=0; j<vPolygons[i].vIndex.size(); ++j)
       {
          out << vPolygons[i].vIndex[j] << " ";
       }
       out << "\n";
    }

//    cout << "vPolygons = " << vPolygons.size() << "\n";

    getline(in,line);  // get rid of the '\n'

    getline(in,line);
//    cout << line << "\n";

	cut_mfix.CalculateCenterOfMass(vPoints,vPolygons);

        bNo_vtk = false;
}

bool CloseTo(double x1 , double y1 , double x2 , double y2)
{
	static const double tolerance = 0.0001;

	if (std::fabs(x1-x2) > tolerance) return false;
	if (std::fabs(y1-y2) > tolerance) return false;

	return true;
}




bool Match2d(double x1 , double y1 , double x2 , double y2 , const CutPolygon & polygon)
{

	double xmin = 1.0e+32;
	double xmax = -xmin;
	double ymin = 1.0e+32;
	double ymax = -ymin;

	for (size_t i=0; i<polygon.vIndex.size(); ++i)
	{
		int k = polygon.vIndex[i];

		if (vPoints[k].x < xmin) xmin = vPoints[k].x;
		if (vPoints[k].x > xmax) xmax = vPoints[k].x;

		if (vPoints[k].y < ymin) ymin = vPoints[k].y;
		if (vPoints[k].y > ymax) ymax = vPoints[k].y;
	}

	double xavg = (xmin+xmax) * 0.5;
	double yavg = (ymin+ymax) * 0.5;

	bool b1 = x1<=xavg && xavg<=x2;	
	bool b2 = y1<=yavg && yavg<=y2;
	
	return b1 && b2;
}

void CutMfix::CalculateCenterOfMass(const vector<CutPoint>   & vPoints ,
		                            const vector<CutPolygon> & vPolygons)
{
	vCenterOfMass.clear();

	for (size_t i=0; i<vPolygons.size(); ++i)
	{
		double x = 0;
		double y = 0;
		double z = 0;

		for (size_t j=0; j<vPolygons[i].vIndex.size(); ++j)
		{
			x +=  vPoints[ vPolygons[i].vIndex[j] ].x;
			y +=  vPoints[ vPolygons[i].vIndex[j] ].y;
			z +=  vPoints[ vPolygons[i].vIndex[j] ].z;
		}

		x /= static_cast<double>( vPolygons[i].vIndex.size() );
		y /= static_cast<double>( vPolygons[i].vIndex.size() );
		z /= static_cast<double>( vPolygons[i].vIndex.size() );

		vCenterOfMass.push_back( CutPoint(x,y,z) );
	}


	const int nx = dx.size();
	const int ny = dy.size();
//	const int nz = dz.size();

//	cout << nx << ":" << ny << ":" << nz << "\n";

        string map_name = Data::sCWD + "cutcell_map.txt";
	ifstream in(map_name.c_str());

	int count = 0;
	if (!in)
	{
	//	cout << "not in\n";
		ofstream mm(map_name.c_str());
		double x1 = -dx[0];
		for (int i=0; i<nx-1; ++i)
		{
	//		x1 += dx[i];
	//		double x2 = x1 + dx[i+1];
			double x2 = x1 + dx[i];
			double y1 = -dy[0];
			for (int j=0; j<ny-1; ++j)
			{
			//	y1 += dy[j];
			//	double y2 = y1 + dy[j+1];
				double y2 = y1 + dy[j];
				bool bFound = false;
				// here ... want to see if override to block exists
				size_t k = 0;
				while (k<vPolygons.size() && !bFound)
				{
					if (Match2d(x1,y1,x2,y2,vPolygons[k]))
					{
						map_IJ_to_vpolygons[ IJ(i,j) ] = k;
						++count;

						set_k.insert(k);
						mm << i << "\t" << j << "\t" << k << "\n";
						bFound = true;
					}
					++k;
				}
				if (!bFound)
				{
				//	mm << i << "\t" << j << "\t" << -1 << "\n";
				}
				y1 += dy[j];
			}
			x1 += dx[i];
		}
	}
	else
	{
//		cout << "in\n";
		int i,j,k;

		while (in >> i >> j >> k)
		{
			map_IJ_to_vpolygons[ IJ(i,j) ] = k;
			set_k.insert(k);
			++count;
		}
	}



//	cout << "count = " << count << "\n";
//	cout <<  map_IJ_to_vpolygons.size() << "\n";
//	cout << "set_k.size() = " << set_k.size() << "\n";

//	set<int>::iterator it = set_k.begin();
//	while (it != set_k.end())
//	{
//		pan << *it << "\n";
//		++it;
//	}

//	{
//		cout << "**************************************\n";
//
//		const CutPolygon & p = vPolygons[0];
//
//		for (int i=0; i<p.vIndex.size(); ++i)
//		{
//			int k = p.vIndex[i];
//
//			cout << vPoints[k].x << "\t" << vPoints[k].y << "\n";
//		}
//
//
//		cout << "**************************************\n";
//	}

}




void set_geom(double * ddx , const int & nx ,
			   double * ddy , const int & ny ,
			   double * ddz , const int & nz )
{
	dx.reserve(nx);
	dy.reserve(ny);
	dz.reserve(nz);

	dx.assign(ddx,ddx+nx);
	dy.assign(ddy,ddy+ny);
	dz.assign(ddz,ddz+nz);

//	cout << nx << ":" << ny << ":" << nz << "\n";

//	for (int i=0; i<ny; ++i)
//	{
//		cout << i << "\t" << dy[i] << "\n";
//	}
}

void plot_ij()
{
	glColor3f(1,0,0);

	for (size_t i=0; i<vPolygons.size(); ++i)
	{
		const CutPolygon & cp = vPolygons[i];

		const vector<int> & v = cp.vIndex;

		glBegin(GL_POLYGON);

		for (size_t k=0; k<v.size(); ++k)
		{
			int j = v[k];

			glVertex3f( vPoints[j].x , vPoints[j].y , 0.0);
		}

		glEnd();
	}
}



  
void plot_ij_slice_(const int & k_use , float * ra , float * ga , float * ba ,
					int * color_array , const int & filling_pass , 
                                        const int * flag , int & status)
{
        if (dz.size() != 1    ||   bNo_vtk)
        {
            status = 0;
            return;
        }

        status = 1;

	const int imax2 = dx.size();
	const int jmax2 = dy.size();
//	pan << "****************************\n";

	static bool bfirst = true;

	int pppp = 0;

	set<int> ssss;

	//int count = (k_use-1) * imax2 * jmax2 - 1;
	int count = (k_use) * imax2 * jmax2 - 1;

	double y1 = -dy[0];
	
	for (int j=0; j<jmax2; ++j)
	{
		double y2 = y1 + dy[j+1];

		double x1 = -dx[0];

		for (int i=0; i<imax2; ++i)
		{
			double x2 = x1 + dx[i+1];
			++count;
			
			if (flag[count] == 1)
			{

			      IJ ij(i,j);

			      map<IJ,int>::iterator it = map_IJ_to_vpolygons.find(ij);


			      if (it == map_IJ_to_vpolygons.end())
			      {			
				      if (filling_pass == 1)
				      {
					      int icol = color_array[count]; // - 1;

					      if (icol != -1)
					      {
						      glColor3f(ra[icol],ga[icol],ba[icol]);
						      ++pppp;
						      glBegin(GL_POLYGON);

						      glVertex3f(x1,y1,0.0f);
						      glVertex3f(x2,y1,0.0f);
						      glVertex3f(x2,y2,0.0f);
						      glVertex3f(x1,y2,0.0f);
						      glVertex3f(x1,y1,0.0f);

						      glEnd();
					      }
				      }
				      
			      }
			      else
			      {
				      int icol = color_array[count]; // - 1;
				      
				      ssss.insert(it->second);

				      if (icol != -1)
				      {
					      if (it->second < (int)vPolygons.size())
					      {
						      glColor3f(1,1,0);
						      glColor3f(ra[icol],ga[icol],ba[icol]);

						      glBegin(GL_POLYGON);
						      CutPolygon cp = vPolygons[it->second];
						      for (size_t i=0; i<cp.vIndex.size(); ++i)
						      {
							      int k = cp.vIndex[i];
							      glVertex3f( vPoints[k].x , vPoints[k].y , 0.0);
						      }
						      glEnd();
					      }
					      else
					      {
					      //	pan << it->second << " : " << vPolygons.size() << "\n";
					      }
				      }
					      
			      }
			}

		

			x1 = x2;
		}

		y1 = y2;
	}

//	pan << "ssss.size() = " << ssss.size() << "\n";
	bfirst = false;
}


  
