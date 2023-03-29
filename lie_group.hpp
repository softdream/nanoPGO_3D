#ifndef __LIE_GROUP_H
#define __LIE_GROUP_H

#include <iostream>
#include <Eigen/Dense>
#include <type_traits>

namespace pgo3d
{

template<typename T>
struct is_numerical_type
{
        static const bool value = false;
};

template<>
struct is_numerical_type<float>
{
        static const bool value = true;
};

template<>
struct is_numerical_type<double>
{
        static const bool value = true;
};

template<>
struct is_numerical_type<uint8_t>
{
	static const bool value = true;
};

template<>
struct is_numerical_type<int8_t>
{
	static const bool value = true;
};

template<>
struct is_numerical_type<uint32_t>
{
        static const bool value = true;
};


template<>
struct is_numerical_type<int32_t>
{
        static const bool value = true;
};

// SO3 Group
template<typename DataType>
class SO3
{
public:
	using Vector3 = typename Eigen::Matrix<DataType, 3, 1>;
	using Matrix3x3 = typename Eigen::Matrix<DataType, 3, 3>;

	template<typename T>
        using Vector3_T = typename Eigen::Matrix<T, 3, 1>;

        template<typename T>
        using Matrix3x3_T = typename Eigen::Matrix<T, 3, 3>;

	SO3()
	{

	}

	SO3( const Matrix3x3& rotaion_matrix ) : rotaion_matrix_( rotaion_matrix )
	{

	}	

	SO3( const Vector3& rotation_vector ) : rotaion_matrix_( rotationVector2Matrix( rotation_vector ) )
	{
		
	}

	template<typename T>
        static 
	const typename std::enable_if_t<is_numerical_type<T>::value, Matrix3x3_T<T>> 
	hat( const Vector3_T<T>& vector ) 
        {
                Matrix3x3_T<T> hat_mat;
                hat_mat << 0, -vector(2), vector(1),
                           vector(2), 0, -vector(0),
                           -vector(1), vector(0), 0;
                return hat_mat;
        }

	template<typename T>
        static
        const typename std::enable_if_t<is_numerical_type<T>::value, Matrix3x3_T<T>>
       	rotationVector2Matrix( const Vector3_T<T>& rotation_vec ) 
        {
                Eigen::AngleAxis<T> rotation_vector( rotation_vec.norm(), rotation_vec.normalized() );

                return rotation_vector.toRotationMatrix();
        }

	template<typename T>
        static
        const typename std::enable_if_t<is_numerical_type<T>::value, SO3<T>>
        exp( const Vector3_T<T>& rotation_vec ) 
        {
		return SO3<T>( rotationVector2Matrix( rotation_vec  ) );
	}


	template<typename T>
	static 
	const typename std::enable_if_t<is_numerical_type<T>::value, Vector3_T<T>>
	rotationMatrix2Vector( const Matrix3x3_T<T>& rotation_mat )
	{
		Eigen::AngleAxis<T> rotaion_vec_eigen;
                rotaion_vec_eigen.fromRotationMatrix( rotation_mat );

                return rotaion_vec_eigen.angle() * rotaion_vec_eigen.axis();

	}

	template<typename T>
	static
	const typename std::enable_if_t<is_numerical_type<T>::value, Vector3_T<T>> 
	log( const SO3<T>& pose_SO3 ) 
	{
		return rotationMatrix2Vector( pose_SO3.getRotationMatrix() );
	}

	const Matrix3x3& getRotationMatrix() const
	{
		return rotaion_matrix_;
	}

	const Matrix3x3 operator * ( const Matrix3x3& other )
	{
		return this->rotaion_matrix_ * other;
	}

private:
	// SO3 Group Element
	Matrix3x3 rotaion_matrix_ = Matrix3x3::Zero();
};

using SO3D = SO3<double>;
using SO3I = SO3<int>;
using SO3F = SO3<float>;

// SE3 Group
template<typename DataType>
class SE3
{
public:
	using Vector3 = typename Eigen::Matrix<DataType, 3, 1>;
	using Vector6 = typename Eigen::Matrix<DataType, 6, 1>;
       	using Matrix3x3 = typename Eigen::Matrix<DataType, 3, 3>;
	using Matrix4x4 = typename Eigen::Matrix<DataType, 4, 4>;
	using Matrix6x6 = typename Eigen::Matrix<DataType, 6, 6>;


	template<typename T>
        using Vector3_T = typename Eigen::Matrix<T, 3, 1>;

	template<typename T>
        using Vector6_T = typename Eigen::Matrix<T, 6, 1>;

	template<typename T>
        using Matrix3x3_T = typename Eigen::Matrix<T, 3, 3>;

        template<typename T>
        using Matrix4x4_T = typename Eigen::Matrix<T, 4, 4>;


	SE3()
	{

	}

	SE3( const Matrix4x4& transform ) : transform_( transform ) 
	{

	}

	SE3( const Vector3& rotation, const Vector3& translation ) : transform_( rotationAndTranslation2Transform( rotation, translation ) )
	{
		std::cout<<"transform_ : "<<std::endl<<transform_<<std::endl;
	}

	// r_t_vec = [ rotation, translation ]
	SE3( const Vector6& r_t_vec ) : transform_( rotationAndTranslation2Transform( r_t_vec ) )	
	{

	}

	template<typename T>
	static
	const typename std::enable_if_t<is_numerical_type<T>::value, Matrix4x4_T<T>> 
	rotationAndTranslation2Transform( const Vector3_T<T>& rotation, const Vector3_T<T>& translation )
	{
		Matrix4x4_T<T> transform = Matrix4x4_T<T>::Zero();
		
		transform.template block<3, 3>( 0, 0 ) = SO3<T>::rotationVector2Matrix( rotation );
		transform.template block<3, 1>( 0, 3 ) = translation;
		transform(3, 3) = 1;

		return transform;
	}

	template<typename T>
	static 
	const typename std::enable_if_t<is_numerical_type<T>::value, Matrix4x4_T<T>> 
	rotationAndTranslation2Transform( const Vector6_T<T>& r_t_vec )
        {
                Matrix4x4_T<T> transform = Matrix4x4_T<T>::Zero();

                transform.template block<3, 3>( 0, 0 ) = SO3<T>::rotationVector2Matrix( r_t_vec.template head<3>() );
                transform.template block<3, 1>( 0, 3 ) = r_t_vec.template tail<3>();
		transform(3, 3) = 1;

                return transform;
        }

	template<typename T>
	static 
	const typename std::enable_if_t<is_numerical_type<T>::value, Matrix3x3_T<T>> 
	caculateJ( const T& theta, const Vector3_T<T>& a )
	{
		auto theta_sin_tmp = ::sin( theta ) / theta;
                auto theta_cos_tmp = ( 1 - ::cos( theta ) ) / theta;

                return theta_sin_tmp * Matrix3x3_T<T>::Identity() + ( 1 - theta_sin_tmp ) * a * a.transpose() + theta_cos_tmp * SO3<T>::hat( a );
	} 

	template<typename T>
        static
        const typename std::enable_if_t<is_numerical_type<T>::value, Matrix3x3_T<T>>
        caculateJInverse( const T& theta, const Vector3_T<T>& a )
        {
		auto half_theta = 0.5 * theta;
		auto theta_cot = ::cos( half_theta ) / ::sin( half_theta );

		return half_theta * theta_cot * Matrix3x3_T<T>::Identity() + ( 1 - half_theta * theta_cot ) * a * a.transpose() - half_theta * SO3<T>::hat( a );
        }


	template<typename T>
	static 
	const typename std::enable_if_t<is_numerical_type<T>::value, Vector6_T<T>> 
	log( const SE3<T>& pose_SE3 ) 
	{
		Matrix3x3_T<T> rotation = pose_SE3.getTransformMatrix().template block<3, 3>( 0, 0 );
		auto translation = pose_SE3.getTransformMatrix().template block<3, 1>( 0, 3 );	

		auto theta = ::acos( ( rotation.trace() - 1 ) * 0.5 );

		auto phi = SO3<T>::rotationMatrix2Vector( rotation );
		auto a = phi.normalized();

		Vector6_T<T> cauthy = Vector6_T<T>::Zero();
		cauthy.template head<3>() = caculateJInverse( theta, a ) * translation;	

		cauthy.template tail<3>() = phi;

		return cauthy;
	}

	template<typename T>
	static 
	const typename std::enable_if_t<is_numerical_type<T>::value, SE3<T>> 
	exp( const Vector6_T<T>& cauthy )
	{
		Vector3_T<T> pho = cauthy.template head<3>(); // J^(-1) * translation
	
		Vector3_T<T> phi = cauthy.template tail<3>(); // rotation
	
		auto theta = phi.norm();
		auto a = phi.normalized();
	
		auto J = caculateJ( theta, a );

		Matrix4x4_T<T> transform = Matrix4x4_T<T>::Zero();
		transform.template block<3, 3>( 0, 0 ) = SO3<T>::rotationVector2Matrix( phi );

		transform.template block<3, 1>( 0, 3 ) = J * pho;
		
		transform(3, 3) = 1;		

		return transform;
	}


	const Matrix4x4& getTransformMatrix() const
        {
                return transform_;
        }

	const Matrix4x4 operator* ( const Matrix4x4& other )
	{
		return this->transform_ * other;
	}
	
	const SE3<DataType> inverse() const
	{
		return SE3<DataType>( this->transform_.inverse() );
	} 

	const Matrix3x3 getRotation() const
	{
		return transform_.template block<3, 3>( 0, 0 );
	}

	const Vector3 getTranslation() const
	{
		return transform_.template block<3, 1>( 0, 3 );
	}

	const Matrix6x6 Adj()
	{
		auto rotation = transform_.template block<3, 3>( 0, 0 );
		auto translation = transform_.template block<3, 1>( 0, 3 );		

		Matrix6x6 ret = Matrix6x6::Zero();
		ret.template block<3, 3>( 0, 0 ) = rotation;
		ret.template block<3, 3>( 3, 3 ) = rotation;
		ret.template block<3, 3>( 0, 3 ) = SO3<DataType>::hat( translation ) * rotation;
		
		return ret;
	}

	const SO3<DataType> getSO3() const
	{
		return SO3<DataType>( transform_.template block<3, 3>( 0, 0 ) );
	}

private:
	// SE3 Group Element
	Matrix4x4 transform_ = Matrix4x4::Zero();
};

using SE3D = SE3<double>;
using SE3F = SE3<float>;
using SE3I = SE3<int>;

}

#endif
