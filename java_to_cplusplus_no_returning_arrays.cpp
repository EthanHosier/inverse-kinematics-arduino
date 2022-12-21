// java_to_cplusplus_no_returning_arrays.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>

/*
  Denavit-Hartenberg Matrix Methods
  https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters

  Each matrix is stored as a 2D vector
*/

//CONSTANTS
const double D_VALS[] = { 0, 85, 290, -20 }; //dimension constants of robot
const double R_VALS[] = { 45, 0, -25, 290 };
const double ALPHAS[] = { M_PI / 2, M_PI / 2, 3 * M_PI / 2, 0 };

//theta angles for aligning x axis: Sigma theta = thetas[n] + anti-clockwise-angle[n]
const double THETAS[] = { 0, 90, 0, -90 };

//used as dx in gradient calculation (is sufficiently small)
const double DX = 0.0001;

//step in degrees per iteration
const double LEARNING_RATE = 0.3;

const int D_H_MATRIX_LENGTH = 4;

void setAsDHMatrix(double array2D[4][D_H_MATRIX_LENGTH] ,double  sigma_theta, double d, double r, double alpha) {

    /*
    Method to set the passed array2D as the corresponding Denavit-Hartenberg matrix
    Matrix produced is the product of the 4 individual transformations:
    M = R(z)*T(z)*T(x)*R(x)

    Theta and alpha received in degrees, thus must be converted to radians
    Both are anti-clockwise relative to reference frame.
    d and r in mm
    */

    //degrees -> radians
    sigma_theta = sigma_theta * M_PI / 180;

    //populating matrix

    //column 0
    array2D[0][0] = cos(sigma_theta);
    array2D[1][0] = sin(sigma_theta);
    array2D[2][0] = 0;
    array2D[3][0] = 0;
    
    //column 1
    array2D[0][1] = -sin(sigma_theta) * cos(alpha);
    array2D[1][1] = cos(sigma_theta) * cos(alpha);
    array2D[2][1] = sin(alpha);
    array2D[3][1] = 0;
    
    //column 2
    array2D[0][2] = sin(sigma_theta) * sin(alpha);
    array2D[1][2] = -cos(sigma_theta) * sin(alpha);
    array2D[2][2] = cos(alpha);
    array2D[3][2] = 0;

    //column 3
    array2D[0][3] = r * cos(sigma_theta);
    array2D[1][3] = r * sin(sigma_theta);
    array2D[2][3] = d;
    array2D[3][3] = 1;
    
}


void setAsProductOf(double product_matrix[][4], double m1[][4], double m2[][4]) {

    /*
    Function which multiplies (dot product) two matrices, setting productMatrix as the result
    */

    double sum;

    //y
    for (int i = 0; i < sizeof(m1); i++) {

        //x
        for (int j = 0; j < sizeof(m1); j++) {

            sum = 0;

            //get dot product sum for each new element
            for (int k = 0; k < sizeof(m1); k++) {
                sum += m1[i][k] * m2[k][j];
            }

            product_matrix[i][j] = sum;
        }
    }

}

void setAsDirectionalJointAngles(double joint_angles[4]) {
    /*
    takes the jointAngles[] , and sets these angles as the anticlockwise angles
    relative to reference frame
    Hard-coded, as each angle being anticlockwise or not already is specific to each joint
    */

    joint_angles[0] = 360.0 - joint_angles[0];

    /*
    These are currently already anticlockwise
    directional_joint_angles[1] = joint_angles[1];
    directional_joint_angles[2] = joint_angles[2];
    directional_joint_angles[3] = joint_angles[3];
    */
    
}


void setAsEndEffectorPosition(double end_effector_position[3],double joint_angles[4]) {
    // Function which sets the vals of the position vector of the end effector in the form[x, y, z]

    //copy each joint_angle to directional_joint_angles, for them to be adjusted relative to anti-clockwise direction
    
    double directional_joint_angles[4];
    for (int i = 0; i < 4; i++)
    {
        directional_joint_angles[i] = joint_angles[i];
    }
    setAsDirectionalJointAngles(directional_joint_angles);

    
    double final_matrix[4][4];
    double temp_matrix_1[4][4];
    double temp_matrix_2[4][4];

    setAsDHMatrix(temp_matrix_1,directional_joint_angles[0] + THETAS[0], D_VALS[0], R_VALS[0], ALPHAS[0]);

    for (int i = 1; i < 4; i++) {
        setAsDHMatrix(temp_matrix_2, directional_joint_angles[i] + THETAS[i], D_VALS[i], R_VALS[i], ALPHAS[i]);
         setAsProductOf(final_matrix,temp_matrix_1, temp_matrix_2);

         //set temp_matrix_1 equal to final_matrix
         for (int j = 0; j < 4; j++)
         {
             for (int k = 0; k < 4; k++)
             {
                 temp_matrix_1[j][k] = final_matrix[j][k];
             }
         }
    }

    //set end_effector_position
    for (int i = 0; i < 3; i++)
    {
        end_effector_position[i] = final_matrix[i][3];
    }
}

double getDistanceBetween(double p1[], double p2[]) {
    /*
    Cost function for gradient descent
    Calculates the magnitude of vector p2->p1
    */

    double distance_between = 0;
    for (int i = 0; i < sizeof(p1); i++) {
        distance_between += pow((p1[i] - p2[i]), 2);

    }

    distance_between = sqrt(distance_between);
    return(distance_between);
}

double getPartialGradient(int n, double current_joint_angles[], double desired_point[]) {
    /*
    Returns the gradient in turns of theta[n] at the current position of the end effector
    (as a result of the current joint angles)
    */

    //setting new_joint_angles
    double new_joint_angles[sizeof(current_joint_angles)];
    for (int i = 0; i < sizeof(current_joint_angles); i++) {
        new_joint_angles[i] = current_joint_angles[i];
    }
    //add dx to angle which partial gradient is in terms of
    new_joint_angles[n] = current_joint_angles[n] + DX;

    double current_point[3];
    setAsEndEffectorPosition(current_point,current_joint_angles);

    double new_point[3];
    setAsEndEffectorPosition(new_point,new_joint_angles);

    //f'[n](x) = "(f(x1,x2,x3 ...xn + dx... xN) - f(x1,x2,x3 ... xn ... xN)) / dx" for dx -> 0
    double gradient = ((getDistanceBetween(new_point, desired_point) - getDistanceBetween(current_point, desired_point)) / DX);

    return gradient;

}


void setAsPartialGradientVector(double partial_gradient_vector[4],double current_joint_angles[4], double desired_point[]) {

    for (int i = 0; i < sizeof(current_joint_angles); i++) {
        partial_gradient_vector[i] = getPartialGradient(i, current_joint_angles, desired_point);
    }
    
}

int main()
{
    /*
   //ONLY TEMPORARY
   double joint_angles[] = {0,0,0,0};

   double end_effector_position[3]; 
   setAsEndEffectorPosition(end_effector_position, joint_angles);

   std::cout << "x: " << end_effector_position[0] << std::endl;
   std::cout << "y: " << end_effector_position[1] << std::endl;
   std::cout << "z: " << end_effector_position[2] << std::endl;
   */

   /*

   //ONLY TEMPORARY
   std::vector<double> joint_angles = {10,20,30,40};

   std::vector<double> end_effector_position = getEndEffectorPosition(joint_angles);

   std::cout << "x: " << end_effector_position[0] << std::endl;
   std::cout << "y: " << end_effector_position[1] << std::endl;
   std::cout << "z: " << end_effector_position[2] << std::endl;

  */
  /*

      joint_angles stores the current angles of each joint IN DEGREES
      Does not account for anti-clockwise or clockwise direction relative to reference frame, thus must be adjusted later
      (0,0,0,0) is at rest, robot in position:

       % ==
      ||
      ||

      NOTE: (0,0,0,0) are not the actual angles of the robot arms - but instead relative to when the reference frames were set.
      Thus, angles are adjusted when actually writing final angles to each servo.
   */

   //angles of the joints at the start [Right_Omoplate, Right_Shoulder, Right_Rotator, Right_Bicep]
    double joint_angles[4] = { 0,0,0,0 }; //THIS IS TO BE CHANGED

    double current_end_effector_position[3];
    setAsEndEffectorPosition(current_end_effector_position,joint_angles);

    //x,y,z (in mm)
    double desired_end_effector_position[] = { -391.5794715386607,197.02194027855472,266.7013179421913 };


    double distance_between = getDistanceBetween(current_end_effector_position, desired_end_effector_position);

    double partial_gradient_vector[sizeof(joint_angles)];

    int iteration = 0;

    while (distance_between > 5) {

        iteration++;

        setAsPartialGradientVector(partial_gradient_vector,joint_angles, desired_end_effector_position);

        for (int i = 0; i < 4; i++) {
            joint_angles[i] -= LEARNING_RATE * partial_gradient_vector[i];
        }

        setAsEndEffectorPosition(current_end_effector_position,joint_angles);

        distance_between = getDistanceBetween(current_end_effector_position, desired_end_effector_position);


        std::cout << "Iteration " << iteration << ":    Distance Between:" << distance_between << std::endl;

    }

    std::cout << std::endl << std::endl;


    for (int i = 0; i < 4; i++) {
        std::cout << "Joint variable " << i << ": " << joint_angles[i] << " degrees" << std::endl;
    }

    std::cout << "\nCurrent Position:" << std::endl;
    std::cout << "x: " << current_end_effector_position[0] << std::endl;
    std::cout << "y: " << current_end_effector_position[1] << std::endl;
    std::cout << "z: " << current_end_effector_position[2] << std::endl;
  

}
