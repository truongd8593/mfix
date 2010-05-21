#ifndef MYFONT_H_ABGTHYJLSJ
#define MYFONT_H_ABGTHYJLSJ

#include <map>
#include <vector>
#include <fstream>

struct MyPoint
{
	int x,y;
};


struct MyFont
{
	std::map<char,std::vector< std::vector<MyPoint> > > char_maps;
	
	float min_x;
	float min_y;
	float max_x;
	float max_y;

	bool FoundCharacter(std::istream & in , char & c);
	void Read(const char * fname);
	void ReadCharacter (std::istream & in, char c);
	void DrawChar(double x , double y , char c) const;
	void DrawString(double x , double y , const char * s) const ;
	void Process();
	void Debug(const char * fname) const;
};




#endif
