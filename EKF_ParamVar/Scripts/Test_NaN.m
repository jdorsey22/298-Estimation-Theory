
thing = isnan(V); 
list  =0 ; 
counter = 1; 

for k =1:length(t)    
    if thing(k) == 1
       
        list(counter) = k;         
        counter = counter +1;        
         
    end    
end 