/* Copyright (C) 2021 Phinite Lab. All rights reserved.
 Contact information
 ------------------------------------------------------
 Developers: 
 Mehmet İşcan    : mehmetmkt@gmail.com
 Ali Burak Özden : ozdenaliburak@gmail.com 
 Company: 
 Phinite Engineering, Istanbul
 Web      :  http://phinitelab.com
 e-mail   :  info@phinitelab.com
 
 */
#ifndef __PHI_MATH_PARAMETERS_HPP__
#define __PHI_MATH_PARAMETERS_HPP__

// in order to generalize our code, we need to use
#ifdef __cplusplus
extern "C"
{
#endif

    ////////////////////////////////////////////////////////////
    // including some libraries for using input/output functions

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

    // including some libraries for using input/output functions
    ////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

////////////////////////////////////////////////////////////
// define error types

#define ALLOCATION_ERROR 1
#define INCONSISTENT_ROW_COLUMN 2
#define FILE_OPEN_ERROR 3
#define SAMPLING_RATE_ERROR 4

// define error types
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// CLASS DEFINITION

class phiMathParameters
{
private:
protected:
public:
    ////////////////////////////////////////////////////////////
    // general variables

    // general variables
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    phiMathParameters()
    {
    }

    ~phiMathParameters()
    {
    }

    // CONSTRUCTORS
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // function prototype decleration

    float **phiMathCreatingEmptyMatrices(int rows, int cols)
    {

        /*
         * creating dynamically empty matrices for general usage
         * output -> address of empty matrices
         * input  -> rows and columns values of matrices
         * */

        float **pd = (float **)malloc(rows * sizeof(float *));

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < rows; i++)
            pd[i] = (float *)malloc(cols * sizeof(float));

        return pd;
    }

    //////////////////////////////////////////////////////////////////////
    // error Handler

    void phiMathPhiErrorHandler(int errorType)
    {
        /*
         * notifying the error result
         * output -> return nothing
         * input  -> type of error
         * */

        switch (errorType)
        {
        case FILE_OPEN_ERROR:
            printf("System Dynamic Parameter files cannot be created!\n");
            break;
        case INCONSISTENT_ROW_COLUMN:
            printf("The rows and columns are not consistent!\n");
            break;
        case ALLOCATION_ERROR:
            printf("Memory allocation cannot be done!\n");
            break;
        case SAMPLING_RATE_ERROR:
            printf("Sampling period cannot be assigned to either negative or zero value!\n");
            break;

        default:
            break;
        }
    }

    // error Handler
    //////////////////////////////////////////////////////////////////////

    float **phiMathPhiVectorMatrixMultiplication(float **firstTerm,
                                                 float **SecondTerm,
                                                 int row1,
                                                 int col1,
                                                 int row2,
                                                 int col2)
    {
        /*
         * creating dynamically matrices produced by multiplication
         * output -> address of multiplied matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of first Matrices
         *           col of first Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        float **pd = this->phiMathCreatingEmptyMatrices(row1, col2);
        float sum = 0;

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row1; i++)
        {
            for (int j = 0; j < col2; j++)
            {
                pd[i][j] = 0;
            }
        }

        for (int i = 0; i < row1; i++) //row of first matrix
        {
            for (int j = 0; j < col2; j++) //column of second matrix
            {
                sum = 0;
                for (int k = 0; k < col1; k++)
                {
                    sum = sum + firstTerm[i][k] * SecondTerm[k][j];
                }
                pd[i][j] = sum;
            }
        }
        return pd;
    }

    float **phiMathPhiSkalarMatrixMultiplication(float skalarTerm,
                                                 float **SecondTerm,
                                                 int row,
                                                 int col)
    {
        /*
         * creating dynamically matrices produced by skalar multiplication
         * output -> address of multiplied matrices
         * input  -> value of skalar term
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        float **pd = this->phiMathCreatingEmptyMatrices(row, col);

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                pd[i][j] = skalarTerm * SecondTerm[i][j];
            }
        }

        return pd;
    }

    float **phiMathPhiMatrixSummation(float **firstTerm,
                                      float **SecondTerm,
                                      int row,
                                      int col)
    {
        /*
         * creating dynamically matrices produced by summation
         * output -> address of summed matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        float **pd = this->phiMathCreatingEmptyMatrices(row, col);

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                pd[i][j] = firstTerm[i][j] + SecondTerm[i][j];
            }
        }

        return pd;
    }

    void phiMathPhiMatrixAssignment(float **assignedTerm,
                                    float **SecondTerm,
                                    int row,
                                    int col)
    {
        /*
         * creating dynamically matrices produced by summation
         * output -> address of summed matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                assignedTerm[i][j] = SecondTerm[i][j];
            }
        }
    }

    // function prototype decleration
    ////////////////////////////////////////////////////////////
};

typedef phiMathParameters *phiMathParametersPtr;

// CLASS DEFINITION
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#endif