#ifndef __NANO_PGO_3D_HPP
#define __NANO_PGO_3D_HPP


#include <Eigen/SparseCholesky>
#include <cmath>

#include <Eigen/Sparse>

#include <vector>
#include <memory>

#include "lie_group.h"

namespace pgo3d
{

template<typename T>
class GraphOptimizer
{
public:
        using DataType = T;
        using Vector3 = typename Eigen::Matrix<DataType, 3, 1>;
        using Vector6 = typename Eigen::Matrix<DataType, 6, 1>;
        using Matrix6x6 = typename Eigen::Matrix<DataType, 6, 6>;
        using Matrix4x4 = typename Eigen::Matrix<DataType, 4, 4>;
        using Matrix3x3 = typename Eigen::Matrix<DataType, 3, 3>;

	using SO3_Type = SO3<DataType>;
	using SE3_Type = SE3<DataType>;

	GraphOptimizer()
	{

	}

	~GraphOptimizer()
	{

	}
	
	void addVertex( const Vector6& pose, const int id ) 
	{
		// 1. pose -> SE3 type
		SE3_Type vertex( pose ); // pose = [ rotation, translation ]
	
		// 2. SE3 type -> se3 type( Lie Group -> Lie Algebra )
		vectex_se3_poses_.push_back( SE3_Type::log( vertex ) );
				
		// 3. vertex id
		vertex_ids_.push_back( id );
	}

	void addEdge( const Vector6& pose_delta,
		      const int from,
		      const int to,
		      const Matrix6x6& info_matrix )   
	{
		SE3_Type delta( pose_delta );
		edge_means_.push_back( SE3_Type::log( delta ) );
		
		edge_from_ids_.push_back( from );
                edge_to_ids_.push_back( to );
		
		info_matrixes_.push_back( info_matrix );
	}

	void execuGraphOptimization( const int max_iterations = 5 )
        {
		// exception check
		if ( vertex_ids_.size() != vectex_se3_poses_.size()  ) {
			std::cerr<<"Vertex Size Error !"<<std::endl;
			return ;
		}		
		if ( ( edge_means_.size() != edge_from_ids_.size() ) || ( edge_from_ids_.size() != edge_to_ids_.size() )  ) {
			std::cerr<<"Edge Size Error !"<<std::endl;
			return ;
		}

		// 1. allocate the memory 
		int vertex_num = vertex_ids_.size();

		allocateMemory( vertex_num );
		
		// 2. Assin the value to the variable 'x_'
		x_ = Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>>( vectex_se3_poses_.data(), vertex_num * 6, 1 );

		// 3. start iteration
		for ( size_t iter = 0; iter < max_iterations; iter ++ ) {
			estimateOnce();
			
		}

		// 4. update
		memcpy( vectex_se3_poses_.data(), x_.data(), vertex_num * 6 * sizeof( DataType ) );
	}		

private:
	void allocateMemory( const int vertex_num )
        {
                H_.resize( vertex_num * 6, vertex_num * 6 );
                b_.resize( vertex_num * 6, 1 );
                delta_x_.resize( vertex_num * 6, 1 );
                x_.resize( vertex_num * 6, 1 );

                std::cout<<"allocate the momory !"<<std::endl;
        }
	
	void estimateOnce()
        {
		// 1. caculate the H & b
		getHessianDerived();	
		
		// 2. 
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<DataType>> solver;
                //Eigen::SparseLLT<Eigen::SparseMatrix<DataType>> solver;
		solver.compute( H_.sparseView() );
		if (solver.info() != Eigen::Success) {
                        std::cerr << "Decomposition Failed !" << std::endl;
                        return;
                }
 
		// 3. using Gaussian Newton Method : delta_x = -H.inverse() * b;
		delta_x_.setZero();	
		delta_x_ = solver.solve( b_ );
		if (solver.info() != Eigen::Success) {
                        std::cerr << "Solving Failed !" << std::endl;
                        return;
                }

		// 4. update 
		x_ += delta_x_;
        }

        void getHessianDerived()
        {
                H_.setZero();
		b_.setZero();

		// for all <e_ij, Omega_ij> do:
                for ( size_t idx = 0; idx < edge_means_.size(); idx ++ ) {
                        int id_i = edge_from_ids_[idx];
                        int id_j = edge_to_ids_[idx];

			// 1. caculate the Jacobians of vertex_i & vertex_j according to the error function	
			getJacobianDerived();
			
			// 2. reset the coefficient parameters
			b_i_.setZero();
                        b_j_.setZero();
                        H_ii_.setZero();
                        H_ij_.setZero();
                        H_ji_.setZero();
                        H_jj_.setZero();

			// 3. recaculate the coefficient parameters
			auto omega = info_matrixes_[idx];
			b_i_ = -J_i_.transpose() * omega * error_;
			b_j_ = -J_j_.transpose() * omega * error_;
			
			H_ii_ = J_i_.transpose() * omega * J_i_;
			H_ij_ = J_i_.transpose() * omega * J_j_;
			//H_ji_ = J_j_.transpose() * omega * J_i_;
			H_jj_ = J_j_.transpose() * omega * J_j_;
			
			// 4. set the Hessian Matrix
			H_.template block<6, 6>( id_i * 6, id_i * 6 ) += H_ii_;
			H_.template block<6, 6>( id_i * 6, id_j * 6 ) += H_ij_;
			H_.template block<6, 6>( id_j * 6, id_j * 6 ) += H_ij_.transpose();
			H_.template block<6, 6>( id_j * 6, id_j * 6 ) += H_jj_;

			// 5. set the b vector
			b_.template block<6, 1>( id_i * 6 ) += b_i_;
			b_.template block<6, 1>( id_j * 6 ) += b_j_;
                }

		// the first vertex is fixed
		H_.template block<6, 6>( 0, 0 ) = Matrix6x6::Identity();
        }
	
	void getJacobianDerived( const int idx, const int id_i, const int id_j )
        {
		error_.setZero();	
		J_i_.setZero();
		J_j_.setZero();		

                auto v_i = vectex_se3_poses_[id_i];
                auto v_j = vectex_se3_poses_[id_j];
                auto z_ij = edge_means_[idx];

		auto vt_i = SE3_Type::exp( v_i );
		auto vt_j = SE3_Type::exp( v_j );
		auto zt_ij = SE3_Type::exp( z_ij );

		// 1. caculate the error
		error_ = SE3_Type::log( zt_ij.inverse() * ( vt_i.inverse() * vt_j ) );
	
		// 2. caculate the Jacobian of vertex_i & vertex_j
        	J_i_ = -getJRInv( vt_i ) + vt_i.inverse().Adj();
		J_j_ = getJRInv( vt_j ) + vt_j.inverse().Adj();
	}

	const Matrix6x6 getJRInv( const SE3_Type& e )
	{
		Matrix6x6 JR = Matrix6x6::Zero();
		
		JR.template block<3, 3>( 0, 0 ) = SO3_Type::hat( SO3_Type::log( e.getSO3() ) );
		JR.template block<3, 3>( 0, 3 ) = SO3_Type::hat( e.getTranslation() );
		//JR.template block<3, 3>( 3, 0 ) = Matrix3x3::Zero();
		JR.template block<3, 3>( 3, 3 ) = JR.template block<3, 3>( 0, 0 );	
		
		return ( JR * 0.5 + Matrix6x6::Identity() );
	}

private:
	std::vector<int> vertex_ids_;
	std::vector<Vector6> vectex_se3_poses_;
        std::vector<int> edge_from_ids_;
        std::vector<int> edge_to_ids_;
        std::vector<Vector6> edge_means_;
        std::vector<Matrix6x6> info_matrixes_;

	Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> H_;
        Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> b_;
        Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> delta_x_;
        Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> x_;

	Matrix6x6 J_i_ = Matrix6x6::Zero();
	Matrix6x6 J_j_ = Matrix6x6::Zero();
	Vector6 error_ = Vector6::Zero();

	Vector6 b_i_ = Vector6::Zero();	
	Vector6 b_j_ = Vector6::Zero();
	Matrix6x6 H_ii_ = Matrix6x6::Zero();
	Matrix6x6 H_ij_ = Matrix6x6::Zero();
	Matrix6x6 H_ji_ = Matrix6x6::Zero();
	Matrix6x6 H_jj_ = Matrix6x6::Zero();
};

}

#endif
