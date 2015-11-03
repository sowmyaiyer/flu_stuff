package multimapper;

public class Gene {
	
	private double snatchReadCounts = 1.0;   //how many snatches map to me?
	private double capReadCounts = 1.0; //how many Cap-seq reads map to me? 
	private double snatchContribution_unnormalized; //how likely am I lose my cap to flu?
	private double snatchContribution_normalized;
	private double snatchPreference;	//how likely is the flu to snatch my cap?
	private double abundance;
	
	public Gene() {
		this(1.0,1.0);
	}

	public Gene(double snatchReadCounts, double capReadCounts) {
		this.snatchReadCounts = snatchReadCounts;
		this.capReadCounts = capReadCounts;
		this.snatchContribution_unnormalized = snatchReadCounts/(snatchReadCounts+capReadCounts);
	}
	
	public void init(double totalSnatchContribution, double totalSnatchReadCounts, double total_cap_plus_snatch_readcounts)
	{
		this.snatchPreference = this.snatchReadCounts/totalSnatchReadCounts;
		this.snatchContribution_normalized = this.snatchContribution_unnormalized/totalSnatchContribution;
		this.abundance = (this.snatchReadCounts + this.capReadCounts)/total_cap_plus_snatch_readcounts;
	}

	public void setSnatchReadCounts(double snatchReadCounts) {
		this.snatchReadCounts = snatchReadCounts;
	}
	
	public double getSnatchReadCounts() {
		return this.snatchReadCounts;
	}

	public final void setCapReadCounts(double capReadCount) {
		this.capReadCounts = capReadCount;
		
	}
	public double getCapReadCounts() {
		return capReadCounts;
	}

	
	public double getAbundance()
	{
		return this.abundance;
	}
	
	public double getSnatchContribution()
	{
		return this.snatchContribution_normalized;
	}

	
	
	public double getSnatchPreference() {
			return this.snatchPreference;
	}


}
