#ifndef __READ_3D_DATA_H
#define __READ_3D_DATA_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iomanip>

#include <Eigen/Dense>


namespace pgo
{

class Read3DData
{
public:

	typedef struct VertexData_
	{
		VertexData_() {}

		int id = -1;
		Eigen::Vector3d pose = Eigen::Vector3d::Zero();
		Eigen::Quaterniond q;
	}VertexData;

	typedef struct EdgeData_
	{
		EdgeData_() {}

		int id_from = -1;
		int id_to = -1;
	
		Eigen::Vector3d pose = Eigen::Vector3d::Zero();
		Eigen::Quaterniond q;
		Eigen::Matrix<double, 21, 1> upper_triangular_vec;
		
	}EdgeData;

	enum DataType 
	{
		vertex_type, 
		edge_type,
		none
	};

	Read3DData()
	{

	}

	~Read3DData()
	{

	}

	bool openFile( const std::string& file_name )
	{
		file.open( file_name.c_str(), std::ifstream::in );
		
		if ( !file.is_open() ) {
			std::cerr<<"Failed to Open the FILE !"<<std::endl;
                        return false;
		}
	
		std::cout<<"Open the FILE !"<<std::endl;
                return true;
	}

	void closeFile()
	{
		return file.close();
	}

	const DataType readOneData( VertexData& vertex, EdgeData& edge  )
	{
		std::string line;
                std::getline( file, line );

                std::istringstream iss( line );

		std::string tag;
		iss >> tag;
	
		// 1. vertex data
		if ( tag.compare( "VERTEX_SE3:QUAT" ) == 0 ) {
			std::string num;
			iss >> num;
			
			// 1.1 id
			vertex.id = std::stoi( num );
		
			// 1.2 pose
			for ( size_t i = 0; i < 3; i ++ ) {
				iss >> num;
				vertex.pose[i] = std::stod( num );
			}

			// 1.3 Quaternion
			Eigen::Vector4d vec;
			for ( size_t i = 0; i < 4; i ++ ) {
				iss >> num;
				vec[i] = std::stod( num );
			}
			vertex.q = Eigen::Quaterniond( vec );
	
			return vertex_type;
		}

		// 2. edge data
		if ( tag.compare( "EDGE_SE3:QUAT" ) == 0 ) {
			std::string num;
                        iss >> num;

			edge.id_from = std::stod( num );

			iss >> num;
			edge.id_to = std::stod( num );
			
			for ( size_t i = 0; i < 3; i ++ ) {
				iss >> num;
				edge.pose[i] = std::stod( num );
			}

			Eigen::Vector4d vec;
                        for ( size_t i = 0; i < 4; i ++ ) {
                                iss >> num;
                                vec[i] = std::stod( num );
                        }
                        edge.q = Eigen::Quaterniond( vec );
		
			for ( size_t i = 0; i < 21; i ++ ) {
                                iss >> num;
                                edge.upper_triangular_vec[i] = std::stod( num );
                        }

			return edge_type;
		}
		
		return none;
	}
	

	const int endOfFile() const
	{
		return file.eof();
	}

	const long getEdgeCount() const
        {
                return edge_count;
        }

        const long getVertexCount() const
        {
                return vertex_count;
        }

private:
	std::ifstream file;

	long edge_count = 0;
        long vertex_count = 0;
};

}

#endif
