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
    // Set/get methods
    std::vector<double> getBetaBins();
    std::vector<double> getBeta();
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
