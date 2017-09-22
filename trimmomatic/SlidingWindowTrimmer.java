package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class SlidingWindowTrimmer extends AbstractSingleRecordTrimmer
{
	private int windowLength;
	private float requiredQuality;
	private float totalRequiredQuality;

	public SlidingWindowTrimmer(String args)
	{
		String arg[]=args.split(":");		
		windowLength=Integer.parseInt(arg[0]);
		requiredQuality=Float.parseFloat(arg[1]);		
		totalRequiredQuality=requiredQuality*windowLength; // Convert to total
	}

    public SlidingWindowTrimmer(int windowLength, float requiredQuality) {
        this.windowLength = windowLength;
        this.requiredQuality = requiredQuality;
        totalRequiredQuality=requiredQuality*windowLength; // Convert to total
    }

    /*
     * @see org.usadellab.trimmomatic.trim.AbstractSingleRecordTrimmer#processRecord(org.usadellab.trimmomatic.fastq.FastqRecord)
     */
	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		int quals[]=in.getQualityAsInteger(true);

		
		if(quals.length<windowLength)
			return null;

		// find the longest high quality fragment
		int total=0;
		for(int i=0;i<windowLength;i++) {
			total+=quals[i];
			//System.out.print(total + ", " + quals[i] + "\n");
		}

		int bestLength, bestStart, bestEnd, fragLength;

		bestLength = 0;
		bestStart = -1;
		bestEnd   = -1;
		fragLength = 0;

		int preBase = 0;	// 0 = low qualty; 1 = high quality

		if(total >= totalRequiredQuality) {
			preBase = 1;
			fragLength++;
		}

		//System.out.print(total + ", " + preBase + ", " + fragLength + "\n");

		int ii;
		for(ii=0; ii<quals.length-windowLength; ii++)
		{
			total=total-quals[ii]+quals[ii+windowLength];
			
			// window start = i + 2;
			// window end   = i + 1 + windowLength
		
			if(total >= totalRequiredQuality)
			{
				fragLength++;
				preBase = 1;	
			} 
			else
			{
				if (fragLength > bestLength)
				{
					bestLength = fragLength;
					bestStart = ii + 1 - fragLength + 1;
					bestEnd   = ii + 1;
				}

				fragLength = 0;
				preBase = 0;
			}
		}

		//System.out.print(bestStart + ", " + bestEnd + ", " + fragLength + "\n");

		if (fragLength > bestLength)
		{
			bestLength = fragLength;
			bestStart = ii + 1 - fragLength + 1;
			bestEnd   = ii + 1;
		}

		//System.out.print(bestStart + ", " + bestEnd + ", " + fragLength + "\n");

		//bestEnd = bestEnd + windowLength - 1;

		if (bestStart == -1 || bestEnd == -1)
			return null;
	
		// extend the best fragment for final start, end ,length
		int finalStart, finalEnd;

		finalStart = bestStart;
		finalEnd = bestEnd + windowLength - 1;

		//System.out.print(finalStart + ", " + finalEnd + ", " + fragLength + "\n");

		// check if the start could be extend
		/*
		if ( bestStart > 1)
		{
			for(int i=bestStart-1; i>1; i--)
			{
				if (quals[i-1] >= requiredQuality)
				{
					finalStart = i;
				}
				else
				{
					break;
				}
			}
		}
		*/

		// check if the start could be shrink
		for(int i=finalStart; i<finalEnd; i++)
		{
			if ( quals[i-1] >= requiredQuality)
			{
				break;
			}
			else
			{
				finalStart=i+1;
			}
		}

		// check if the end could be extend
		/*
		if ( bestEnd < quals.length )
		{
			for(int i=bestEnd; i<quals.length; i++)
			{
				if (quals[i] >= requiredQuality)
				{
					finalEnd = i+1;
				}
				else
				{
					break;
				}
			}
		}
		*/

		// check if the the end could be shrink
		for(int i=finalEnd; i>finalStart; i--)
		{
			if (quals[i-1] >= requiredQuality)
			{
				break;
			}
			else
			{
				finalEnd = i-1;
			}
		}

		int finalLength = finalEnd - finalStart + 1;

		//System.out.print(finalStart + ", " + finalEnd + ", " + fragLength + "\n");

		// return the trimmed base
		if(finalLength < 1)
			return null;
		
		if(finalLength < quals.length)
			return new FastqRecord(in,finalStart-1,finalLength);
		
		return in;
	}
}
