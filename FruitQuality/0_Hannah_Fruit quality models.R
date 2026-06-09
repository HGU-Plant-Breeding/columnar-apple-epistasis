

######JA
# try first:
trait ~ pop + rootstock,
Random=  ~ diag(loc_yr):vm(Genotype, Ginv) + at(loc_yr):row + at(loc_yr):column + year,  
Residual: ~ dsum(~ ar1(row):ar1(column) | loc_yr)                                        

# try second:

trait ~ pop + rootstock,
Random=  ~ corh(loc_yr):vm(Genotype, Ginv) + at(loc_yr):row + at(loc_yr):column + year,  
Residual: ~ dsum(~ ar1(row):ar1(column) | loc_yr)  

# try third:

trait ~ pop + rootstock,
Random=  ~ us(loc_yr):vm(Genotype, Ginv) + at(loc_yr):row + at(loc_yr):column + year,  
Residual: ~ dsum(~ ar1(row):ar1(column) | loc_yr) 


##Back up
trait ~ pop + rootstock,
Random=  ~ fa(loc_yr, 1):vm(Genotype, Ginv) + at(loc_yr):row + at(loc_yr):column + year,  
Residual: ~ dsum(~ ar1(row):ar1(column) | loc_yr)



#####Firmness+RGB

Fixed:    trait ~ pop + rootstock,
Random:   ~ diag(loc_yr):vm(Genotype, Ginv) + diag(loc_yr):row + diag(loc_yr):column + year,
Residual: ~ dsum(~ idv(units)| loc_yr)    

###continue with variance structure modelling 