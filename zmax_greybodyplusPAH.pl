use lib "$ENV{HOME}/Software";
$ENV{DATADIR} = "$ENV{HOME}/Software/filters";

use PDL;
use PDL::NiceSlice;
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;

$ABflag = 0; #use Vega K mag to compare for zmax
$OmegaM = 0.3; $H0 = 70;
$Law = "SMC"; # Dust Law MW/SMC
$nebfac = 2.0; # Extra factor of extinction for neb lines
$pc = 3.086E16; # One parsec in meters
$flux0 = 3631; # AB zero-point in Jy
$clight = 2.99792e8; # speed of light in  meters
$planck = 6.6260755e-34; #planck's constant in m^2 kg /s
$boltzmann = 1.3806503e-23; # m^2 kg s-2 K


open IN, "$ENV{HOME}/Software/Masses/GDDSVIzK1234_greybody_free/SEDparameters.dat" or die "Can't open SED file";
open OUT, ">$ENV{HOME}/Software/Masses/GDDSVIzK1234_greybody_free/zmax.dat" or die "can't open zmax file"; 
print OUT "# id      zmax\n";
select OUT; $| = 1; select STDOUT; $| = 1;

while (<IN>) {
    next if /^#/;
    @v = split;
    $id = $v[0]; $zsp = $v[1]; $specfile = $v[2]; $itime = $v[3]; $fac = $v[4]; $facburst = $v[5]; $AVV =$v[6]; $bbtemp = $v[7];

    print "Getting SED from $specfile\n";

    ($t2,$dummy, $wav,$spec,$emspec) = read_peg_spec($specfile,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A
    $spec *= 1E10; $emspec *= 1E10; # Convert to 1E10 mass galaxy to avoid rounding
         
    $spec2 = $spec->dice_axis(1,$itime); 
    $emspec2 = $emspec->dice_axis(1,$itime); 
           
    $attn = peidust($Law, $wav)/peidust($Law, 5500);
    $spec2 = $spec2 * 10**( -0.4 * $attn * $AVV ) + $emspec2 * 10**(-0.4 * $attn * $nebfac * $AVV);

    
#now fit the additional starburst component. Must interpolate burst spectrum  to match PEGASE SED
    $wav2 = 1e-9 * sequence(20000) + 1e-7; #wav in meters 
    $specburst2 = 2 * $planck * ($clight ** 2) * $wav2 ** (-6) * ( exp( $clight * $planck /($wav2 * $boltzmann * $bbtemp))) ** (-1); #in W/m^3
    
    $specburst2 *= 1e10; #in W/m^2/A 
    
    $wav2 *= 1e10; #convert to angstroms
    
# add 3.3 um PAH feature following Lorentian Profile parametrized by Verstraete et al. (2001)
    
#    $PAH = 2.27 * $specburst2->index(which( $wav2 == 40000)) / (1 + ((1 / ($wav2 * 1e-8) - 3039.1)/19.4)**2); 
    
#    $specburst2 += $PAH;
     
       
    $specburst = interpol($wav,$wav2,$specburst2);     

    $spectotal = $fac*$spec2 + $facburst*$specburst;

    $z = 0;
    $K = 0;

    while ($K <= 20.6) {
	$z+= 0.001;
	($wz, $sz) = redshift_spectra($wav, $spectotal, pdl($z));
	
	# Coerce 1D             
	die "Something wrong \n" if (nelem($sz) ne $sz->getdim(0)) or (nelem($wz) ne $wz->getdim(0)); # Should never happen! 
	$wz = $wz->clump(-1)->copy; $sz = $sz->clump(-1)->copy;
	
	$DL = lumdist($z);
	$sz /= 4*$Pi * ($DL*1E6*$pc)**2; 
	
	# $sz now matches that used to produce the colors for this model
	
	$K = mag($wz,$sz,'sys_circi_K');
	
    }
    $zmax = $z;
    $zmax = 0.0 if ($zmax < $zsp);
    printf OUT "%10s %5.2f\n", $id,$zmax;
}
close IN; close OUT;
