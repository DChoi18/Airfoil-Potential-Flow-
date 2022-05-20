function w = KarmanTrefftzTransform(z,b,k)
w = k*b*((z+b).^k+(z-b).^k)./((z+b).^k-(z-b).^k); %Karman Trefftz Transformation
end