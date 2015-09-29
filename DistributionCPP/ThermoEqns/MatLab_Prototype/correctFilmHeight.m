function scalars = correctFilmHeight(scalars) 
% Function to correct non-physicalities in film height

X = scalars.X_;
Z = scalars.Z_;
mimp = scalars.mimp_;

tolIMAG = 1e-6;
tolNEG = (5e-3)*max(X);
indX = find((abs(imag(X)) > tolIMAG) | (isnan(X)) | (isinf(X)) | X<=tolNEG);
if (~isempty(indX))
    % Set ice accretion rate (z) equal to mimp everywhere where hf = 0
    indFix = indX;
    Z(indFix) = mimp(indFix);
    scalars.Z_ = Z;
    % Reset new physically consistent film height (= 0)
    X(indFix) = 0;
    scalars.X_ = X;
end


end