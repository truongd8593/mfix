#include <QtOpenGL>

#include "myfont.h"
#include "some_font.h"

#include <sstream>

bool MyFont::FoundCharacter(std::istream & in , char & c)
{
	using namespace std;

	string line;
	int    ic;

	while ( getline(in,line) )
	{
		string::size_type pos = line.find("static");

		if (pos != string::npos)
		{
			pos = line.find("F0S");
			if (pos != string::npos)
			{
				stringstream ss(line.substr(pos+3));
				ss >> ic;
				c = static_cast<char>(ic);
				return true;
			}
		}
	}
	return false;
}

void MyFont::ReadCharacter (std::istream & in, char c)
{
	using namespace std;

	string line;

	vector< vector<MyPoint> > vv;
	vector<MyPoint>           v;
	MyPoint                   p;
	char                      comma;

	while (getline(in,line))
	{
		if (line.find("GL_LINE_STRIP") != string::npos)
		{
			getline(in,line);  // ,
			getline(in,line);
			const int N = Convert<string,int>(line);

			for (int i=0; i<N; ++i)
			{
				getline(in,line);  // ,
				getline(in,line);
				stringstream ss(line);

				ss >> p.x >> comma >> p.y;
				v.push_back(p);
				
				if (p.x < min_x) min_x = p.x;
				if (p.x > max_x) max_x = p.x;
				if (p.y < min_y) min_y = p.y;
				if (p.y > max_y) max_y = p.y;
			}
			vv.push_back(v);
			v.clear();
		}
		else if (line.find("};") != string::npos)
		{
			char_maps[c] = vv;
			return;
		}
	}
}

void MyFont::Read(const char * fname)
{
	min_x = 100000000;
	min_y = 100000000;
	max_x = -min_x;
	max_y = -min_y;
	
	std::ifstream in(fname);

	char c;

	while (FoundCharacter(in,c))
	{
		ReadCharacter(in,c);
	}
	
//		std::cout << min_x << "\n";
//		std::cout << min_y << "\n";
//		std::cout << max_x << "\n";
//		std::cout << max_y << "\n";
}

void MyFont::DrawChar(double x , double y , char c) const 
{
	using namespace std;

	map<char,vector< vector<MyPoint> > >::const_iterator it = char_maps.find(c);

	if (it == char_maps.end()) return;

	const vector< vector<MyPoint> > & vv = it->second;

	for (size_t i=0; i<vv.size(); ++i)
	{
		glBegin(GL_LINE_STRIP);
		  for (size_t j=0; j<vv[i].size(); ++j)
		  {
			  glVertex3f( x+vv[i][j].x , y+vv[i][j].y , 0);
		  }
		glEnd();
	}
}

void MyFont::DrawString(double x , double y , const char * s) const 
{
	using namespace std;

	const int N = strlen(s);

	for (int k=0; k<N; ++k)
	{
		DrawChar(x,y,s[k]);
		x += 30000;
	}
}

void MyFont::Process()
{
  using namespace std;
  
//	  ofstream out("process_out.txt");
  
  
  vector< vector<MyPoint> > vv;
  vector<MyPoint>           v;
  MyPoint                   p;
  int                      c;
  int                       n , nv;
  
  SomeFont sf;
  string s = sf.init();
//	  out << s << "\n";
//	  out.close();
  stringstream ss(s);
  
  while (ss >> c >> n)
  {
    for (int i=0; i<n; ++i)
    {
      ss >> nv;
      for (int j=0; j<nv; ++j)
      {
	ss >> p.x >> p.y;
	v.push_back(p);
      }
      vv.push_back(v);
      v.clear();
    }
    char_maps[(char)c] = vv;
    vv.clear();
  }
  
  
}

  

void MyFont::Debug(const char * fname) const
{
	using namespace std;

	ofstream out(fname);
	
	out << "struct SomeFont\n";
	out << "{\n";
	out << "void init()\n";
	out << "{\n";
	out << "\tstringstream ss;\n";

	map<char,vector< vector<MyPoint> > >::const_iterator it = char_maps.begin();

	while (it != char_maps.end())
	{
		out << "\tss << \"" << (int)it->first << " ";
		const vector< vector<MyPoint> > & vv = it->second;
		out << vv.size() << "\" << endl;\n";

		for (size_t i=0; i<vv.size(); ++i)
		{
			out << "\t\tss << " << vv[i].size() << " << endl;\n";
			for (size_t j=0; j<vv[i].size(); ++j)
			{
			  out << "\t\t\tss << \"" << vv[i][j].x << "   " <<  vv[i][j].y << "\" << endl;\n";
			}
		}
		++it;
	}
	out << "}\n";
	out << "};\n";
}




