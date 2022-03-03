function out = getPlanetName(ID)
%GETPLANETNAME Returns Char of Planet Given Numerical ID

    p = {'Mercury';
         'Venus';
         'Earth';
         'Mars';
         'Jupiter';
         'Saturn';
         'Uranus';
         'Neptune';
         'Pluto'};
    out = p{ID};
    
end