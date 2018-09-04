double cech_scale()
{
    //steps for the algorithm
    //read and validate disk system from text file
    //calculate vietori rips scale of the disk system
    //verify if rho is positive (if the minimum output is positive)
    //  if so, then output the vietori rips scale and point of intersection and end
    //if is not completely positive, apply bisection to a 3-disk system rho to approximate cech scale. Repeat for all possible 3-disk systems
    //take the max cech scale calculated
    //
    //validation: if 2d then apply all above steps with any number of disks.
    //            if Rd, R>2 then must be only 3 disks and must apply a projection plane using the 3 centers of the disks, apply above steps and redimention the intersection points to Rd.
    //
    //output cech scale, point of intersection and wheter cech and vieti rips scale coincide
}
