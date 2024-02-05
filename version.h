#ifndef VERSION_H
#define VERSION_H

namespace AutoVersion{
	
	//Date Version Types
	static const char DATE[] = "22";
	static const char MONTH[] = "11";
	static const char YEAR[] = "2019";
	static const char UBUNTU_VERSION_STYLE[] =  "19.11";
	
	//Software Status
	static const char STATUS[] =  "Alpha";
	static const char STATUS_SHORT[] =  "a";
	
	//Standard Version Type
	static const long MAJOR  = 1;
	static const long MINOR  = 0;
	static const long BUILD  = 60;
	static const long REVISION  = 370;
	
	//Miscellaneous Version Types
	static const long BUILDS_COUNT  = 80;
	#define RC_FILEVERSION 1,0,60,370
	#define RC_FILEVERSION_STRING "1, 0, 60, 370\0"
	static const char FULLVERSION_STRING [] = "1.0.60.370";
	
	//These values are to keep track of your versioning state, don't modify them.
	static const long BUILD_HISTORY  = 60;
	

}
#endif //VERSION_H
