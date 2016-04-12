#ifndef __AIRFOIL_H__
#define __AIRFOIL_H__

#include <eigen3/Eigen/Dense>
#include <QuadTree/Bucket.h>
#include <gsl/gsl_histogram.h>
#include <Grid/PLOT3D.h>

class Airfoil {
  public:
    Airfoil(std::vector<double>& X, std::vector<double>& Y);
    ~Airfoil();
    void findPanel(std::vector<double>& XYq, std::vector<double>& XYnn, std::vector<double>& NxNy, std::vector<double>& TxTy);
    void findPanel(std::vector<double>& XYq, std::vector<double>& XYnn, std::vector<double>& NxNy, std::vector<double>& TxTy, int& indexNN);
    double calcIncidenceAngle(std::vector<double>& XYq,std::vector<double>& UVq,int indNN);
    double interpXYtoS(std::vector<double>& XYq);
    void appendFilm(double sCoord, double mass);
    void calcCollectionEfficiency(double fluxFreeStream,double dS);
    void calcStagnationPt(PLOT3D& grid);
    // Methods for updating grid based on thermodynamic ice calculation
    double computeJaggednessCriterion(int id1, int i2, int id3, int id4);
    void correctJagged(int id1, int id2, int id3, int id4);
    std::vector<double> tridiagSolve(std::vector<std::vector<double>>& A, std::vector<double>& r);
    std::vector<std::vector<double>> LaplacianMatrix(int N, double eps);
    std::vector<double> movingAverage(std::vector<double>& X, double smooth);
    void growIce(std::vector<double>& sTHERMO, std::vector<double>& mice, double DT, double chord, const char* strSurf);
    // Set/get methods
    std::vector<double> getBetaBins();
    std::vector<double> getBeta();
    std::vector<double> getX();
    std::vector<double> getY();
    void setStagPt(double sLoc);
    double getStagPt();

  private:
    Eigen::VectorXd panelX_;
    Eigen::VectorXd panelY_;
    Eigen::VectorXd panelS_;
    Eigen::MatrixXd tangent_;
    Eigen::MatrixXd normal_;
    std::vector<double> FilmScoords_;
    std::vector<double> FilmMass_;
    std::vector<double> BetaBins_;
    std::vector<double> Beta_;
    double impingeLimitUP_;
    double impingeLimitDOWN_;
    double stagPt_;
    double stagPtX_;
    double stagPtY_;
    Bucket panelSearcher_;
    void calcSCoords();
    
};



#endif
