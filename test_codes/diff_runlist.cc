#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <sstream>
#include <iostream>
#include <limits>

using namespace std;

int main()
{

   const char * crosslist = "run_list_both.tsv";
   std::cout <<crosslist <<  endl;
   ifstream  fileP(crosslist);


      
   string label;
   getline(fileP, label, '\n');
      
   istringstream gccNeedsThisOnASeperateLine (label);
   const vector<string> Labels( (istream_iterator<string>( gccNeedsThisOnASeperateLine )), (istream_iterator<string>()) );
   
      //Read the remainer and parse it to a 2d vector of doubles
   std::vector <std::vector<double> > Parameters;
   do {
      vector<double> input(Labels.size());
      for(int i = 0; i < input.size(); i++){
         fileP >> input[i];
      }
      if(!fileP.fail()) //control for empty line at end of file
         Parameters.push_back(input);
   } while(fileP.ignore(numeric_limits<streamsize>::max(), '\n'));

     
   std::cout<< " The size of the label/column "<<Labels.size()<<"\n";
   
   for (int i=0;i<Labels.size();i++){
      cout << Labels[i] << "\t" << Parameters[0][i]  <<"\n";
   }
  
   
   for (int i=0;i<Labels.size();i++){
      cout << Labels[i] << "\t" << Parameters[5][i]  <<"\n";
   }
  
   
   int sssize = Parameters.size();
   cout << "The sssize is" << sssize << endl;
   for (int i=0;i<Labels.size();i++){
      cout << Labels[i] << "\t" << Parameters[sssize-1][i]  <<"\n";
   }
   
   int row = 0;
   std::cout<<" what is the value of this thing "<<Parameters[row].size()<<"\n";
   for(int i = 0; i<Parameters[row].size(); i++)
      {
      for(int s = 0; s<Parameters.size(); s++)
         {
          std::cout<<"For i= "<<i<< "and s= "<<s<< " print parameter[s][i]= "<<Parameters[s][i]<<"\n";
         }
      row++;
      }
   
   
   
   return 0;
}
	

