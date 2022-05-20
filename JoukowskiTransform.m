function w = JoukowskiTransform(z,b,direction)
switch direction
    case 'forward'
        w = z+b^2./(z); %Joukowski Transformation
    case 'inverse'
        w = 0.5*(z+sqrt((z).^2-4*b^2));
end

end