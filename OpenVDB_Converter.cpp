// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <vector>

// using namespace std;
// //delete space, tab, etc.
// string Trim(string& str)
// {
// 	str.erase(0, str.find_first_not_of(" \t\r\n"));
// 	str.erase(str.find_last_not_of(" \t\r\n") + 1);
// 	return str;
// }

// int main()
// {
    
// 	ifstream fin("data.csv"); //open
// 	string line;
// 	getline(fin, line);
// 	while (getline(fin, line))   
// 	{
// 		cout << "original:\t" << line << endl; 
// 		istringstream sin(line); 
// 		vector<string> fields; 
// 		string field;
// 		while (getline(sin, field, ',')) 
// 		{
// 			fields.push_back(field); 
// 		}
// 		string distance = Trim(fields[0]); 
// 		string x = Trim(fields[1]); 
// 		string y = Trim(fields[2]); 
// 		string z = Trim(fields[3]);
// 		cout << "processed:\t" << distance << x << "\t" << y << "\t" << z << endl;
// 	}
// 	return EXIT_SUCCESS;
// }

#include <stdio.h>
#include <dirent.h>

int main(int argc, const char**argv) {
    struct dirent *entry = nullptr;
    DIR *dp = nullptr;

    dp = opendir(argc > 1 ? argv[1] : "/");
    if (dp != nullptr) {
        while ((entry = readdir(dp)))
            printf ("%s\n", entry->d_name);
    }

    closedir(dp);
    return 0;
}


// openvdb::FloatGrid::Ptr mygrid = openvdb::FloatGrid::create();
// openvdb::FloatGrid::Accessor accessor = mygrid->getAccessor();
// grid.iterateRegionSerial(cell_region, [&](const IV& cell, TopoGridState<T, dim>& cell_state) {
//     if (!cell_state.inside_cell) return;
//     for (int x = 0; x < 2; ++x)
//         for (int y = 0; y < 2; ++y)
//             for (int z = 0; z < 2; ++z) {
//                 if (remove_empty && rho(cell_state.cell_idx * ppc + 4 * x + 2 * y + z) == 0.) continue;
//                 openvdb::Coord xyz(cell[0] * 2 + x, cell[1] * 2 + y, cell[2] * 2 + z);
//                 accessor.setValue(xyz, float(rho(cell_state.cell_idx * ppc + 4 * x + 2 * y + z)));}});
// // Set the name of the grid
// mygrid->setName("density");
// // Create a VDB file object and write out the grid.
// openvdb::io::File(filename).write({ mygrid });