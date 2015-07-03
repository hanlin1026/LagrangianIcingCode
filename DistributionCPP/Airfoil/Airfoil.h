#ifndef __AIRFOIL_H__
#define __AIRFOIL_H__

#include <eigen3/Eigen/Dense>
#include <QuadTree/Bucket.h>

class Airfoil {
  public:
    Airfoil(Eigen::VectorXd& X, Eigen::VectorXd& Y);
    ~Airfoil();
    void findPanel(std::vector<double>& XYq, std::vector<double>& XYnn, std::vector<double>& NxNy, std::vector<double>& TxTy);
    void findPanel(std::vector<double>& XYq, std::vector<double>& XYnn, std::vector<double>& NxNy, std::vector<double>& TxTy, int& indexNN);
    double calcIncidenceAngle(std::vector<double>& XYq,std::vector<double>& UVq,int indNN);
    double interpXYtoS(std::vector<double>& XYq);
    void appendFilm(double sCoord, double mass);

  private:
    Eigen::VectorXd panelX_;
    Eigen::VectorXd panelY_;
    Eigen::VectorXd panelS_;
    Eigen::MatrixXd tangent_;
    Eigen::MatrixXd normal_;
    std::vector<double> FilmScoords_;
    std::vector<double> FilmMass_;
    double impingeLimitUP_;
    double impingeLimitDOWN_;
    double stagPt_;
    Bucket panelSearcher_;
    void calcSCoords();
    
};



#endif
