#!/bin/bash

# Convert Bismark cx_report to BigWig for IGV visualization

if [ $# -eq 0 ]; then
    echo "Usage: $0 [--at_rename] <cx_report_file> [output_prefix] [chromosome]"
    echo "Example: $0 sample.CX_report.txt sample"
    echo "Example with chromosome filter: $0 sample.CX_report.txt sample Chr1"
    echo "Example with AT chromosome renaming: $0 --at_rename sample.CX_report.txt sample"
    echo "Available chromosomes: Chr1, Chr2, Chr3, Chr4, Chr5, ChrC, ChrM"
    echo ""
    echo "Options:"
    echo "  --at_rename  Rename NCBI RefSeq chromosome names to standard AT names"
    echo "               NC_003070.9 -> Chr1, NC_003071.7 -> Chr2, etc."
    exit 1
fi

# Parse command line arguments
AT_RENAME=false
if [ "$1" = "--at_rename" ]; then
    AT_RENAME=true
    shift
fi

CX_FILE=$1
PREFIX=${2:-$(basename $CX_FILE .CX_report.txt)}
CHROMOSOME=$3

# Check if input file exists
if [ ! -f "$CX_FILE" ]; then
    echo "Error: Input file $CX_FILE not found"
    exit 1
fi

echo "Processing $CX_FILE..."

# Create chromosome renaming function if AT_RENAME is enabled
if [ "$AT_RENAME" = true ]; then
    echo "AT chromosome renaming enabled"
    RENAME_AWK='
    {
        if ($1 == "NC_003070.9") $1 = "Chr1"
        else if ($1 == "NC_003071.7") $1 = "Chr2"
        else if ($1 == "NC_003074.8") $1 = "Chr3"
        else if ($1 == "NC_003075.7") $1 = "Chr4"
        else if ($1 == "NC_003076.8") $1 = "Chr5"
        else if ($1 == "NC_037304.1") $1 = "ChrM"
        else if ($1 == "NC_000932.1") $1 = "ChrC"
    }'
else
    RENAME_AWK=''
fi

# Add chromosome filter message if specified
if [ ! -z "$CHROMOSOME" ]; then
    echo "Filtering for chromosome: $CHROMOSOME"
    PREFIX="${PREFIX}_${CHROMOSOME}"
    
    # Debug: Check if chromosome exists in the file before processing
    echo "Checking for chromosome $CHROMOSOME in file..."
    if [ "$AT_RENAME" = true ]; then
        # When AT_RENAME is enabled, check for both original NCBI names and renamed versions
        chr_count=$(awk -v chr="$CHROMOSOME" '{
            orig_chr = $1
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
            if ($1==chr || orig_chr==chr) count++
        } END {print count+0}' $CX_FILE)
    else
        chr_count=$(awk -v chr="$CHROMOSOME" '$1==chr {count++} END {print count+0}' $CX_FILE)
    fi
    echo "Found $chr_count lines for chromosome $CHROMOSOME"
    
    if [ $chr_count -eq 0 ]; then
        echo "Warning: No data found for chromosome $CHROMOSOME"
        echo "Available chromosomes in file:"
        if [ "$AT_RENAME" = true ]; then
            awk '{
                orig_chr = $1
                if ($1 == "NC_003070.9") $1 = "Chr1"
                else if ($1 == "NC_003071.7") $1 = "Chr2"
                else if ($1 == "NC_003074.8") $1 = "Chr3"
                else if ($1 == "NC_003075.7") $1 = "Chr4"
                else if ($1 == "NC_003076.8") $1 = "Chr5"
                else if ($1 == "NC_037304.1") $1 = "ChrM"
                else if ($1 == "NC_000932.1") $1 = "ChrC"
                print orig_chr " -> " $1
            }' $CX_FILE | sort | uniq -c | head -10
        else
            awk '{print $1}' $CX_FILE | sort | uniq -c | head -10
        fi
        exit 1
    fi
fi

# Create bedGraph files for each methylation context
if [ ! -z "$CHROMOSOME" ]; then
    echo "Processing CG methylation for $CHROMOSOME..."
    if [ "$AT_RENAME" = true ]; then
        awk -v chr="$CHROMOSOME" 'BEGIN{OFS="\t"} {
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
        } $1==chr && $6=="CG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CG.bedGraph
    else
        awk -v chr="$CHROMOSOME" 'BEGIN{OFS="\t"} $1==chr && $6=="CG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CG.bedGraph
    fi
    
    echo "Processing CHG methylation for $CHROMOSOME..."
    if [ "$AT_RENAME" = true ]; then
        awk -v chr="$CHROMOSOME" 'BEGIN{OFS="\t"} {
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
        } $1==chr && $6=="CHG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHG.bedGraph
    else
        awk -v chr="$CHROMOSOME" 'BEGIN{OFS="\t"} $1==chr && $6=="CHG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHG.bedGraph
    fi
    
    echo "Processing CHH methylation for $CHROMOSOME..."
    if [ "$AT_RENAME" = true ]; then
        awk -v chr="$CHROMOSOME" 'BEGIN{OFS="\t"} {
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
        } $1==chr && $6=="CHH" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHH.bedGraph
    else
        awk -v chr="$CHROMOSOME" 'BEGIN{OFS="\t"} $1==chr && $6=="CHH" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHH.bedGraph
    fi
else
    if [ "$AT_RENAME" = true ]; then
        awk 'BEGIN{OFS="\t"} {
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
        } $6=="CG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CG.bedGraph
        
        awk 'BEGIN{OFS="\t"} {
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
        } $6=="CHG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHG.bedGraph
        
        awk 'BEGIN{OFS="\t"} {
            if ($1 == "NC_003070.9") $1 = "Chr1"
            else if ($1 == "NC_003071.7") $1 = "Chr2"
            else if ($1 == "NC_003074.8") $1 = "Chr3"
            else if ($1 == "NC_003075.7") $1 = "Chr4"
            else if ($1 == "NC_003076.8") $1 = "Chr5"
            else if ($1 == "NC_037304.1") $1 = "ChrM"
            else if ($1 == "NC_000932.1") $1 = "ChrC"
        } $6=="CHH" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHH.bedGraph
    else
        awk 'BEGIN{OFS="\t"} $6=="CG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CG.bedGraph
        awk 'BEGIN{OFS="\t"} $6=="CHG" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHG.bedGraph
        awk 'BEGIN{OFS="\t"} $6=="CHH" && $4+$5>0 {print $1,$2-1,$2,$4/($4+$5)}' $CX_FILE | sort -k1,1 -k2,2n > ${PREFIX}_CHH.bedGraph
    fi
fi

# Check if bedGraph files were created successfully
for context in CG CHG CHH; do
    if [ -f ${PREFIX}_${context}.bedGraph ]; then
        lines=$(wc -l < ${PREFIX}_${context}.bedGraph)
        echo "Created ${PREFIX}_${context}.bedGraph with $lines lines"
        if [ $lines -eq 0 ]; then
            echo "Warning: ${context} bedGraph file is empty!"
        fi
    else
        echo "Error: Failed to create ${PREFIX}_${context}.bedGraph"
    fi
done

# Create chromosome sizes file (for Arabidopsis TAIR10)
if [ ! -z "$CHROMOSOME" ]; then
    # Create chromosome sizes file for specific chromosome only
    case "$CHROMOSOME" in
        "Chr1") echo -e "Chr1\t30427671" > ${PREFIX}_chrom.sizes ;;
        "Chr2") echo -e "Chr2\t19698289" > ${PREFIX}_chrom.sizes ;;
        "Chr3") echo -e "Chr3\t23459830" > ${PREFIX}_chrom.sizes ;;
        "Chr4") echo -e "Chr4\t18585056" > ${PREFIX}_chrom.sizes ;;
        "Chr5") echo -e "Chr5\t26975502" > ${PREFIX}_chrom.sizes ;;
        "ChrC") echo -e "ChrC\t154478" > ${PREFIX}_chrom.sizes ;;
        "ChrM") echo -e "ChrM\t366924" > ${PREFIX}_chrom.sizes ;;
        *) echo "Error: Unknown chromosome $CHROMOSOME"; exit 1 ;;
    esac
else
    # Create full chromosome sizes file
    cat > ${PREFIX}_chrom.sizes << EOF
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
ChrC	154478
ChrM	366924
EOF
fi

# Convert bedGraph to BigWig
for context in CG CHG CHH; do
    if [ -s ${PREFIX}_${context}.bedGraph ]; then
        echo "Converting ${context} to BigWig..."
        bedGraphToBigWig ${PREFIX}_${context}.bedGraph ${PREFIX}_chrom.sizes ${PREFIX}_${context}.bw
        rm ${PREFIX}_${context}.bedGraph
    fi
done

rm ${PREFIX}_chrom.sizes

if [ ! -z "$CHROMOSOME" ]; then
    echo "Conversion complete for chromosome $CHROMOSOME. Output files: ${PREFIX}_CG.bw, ${PREFIX}_CHG.bw, ${PREFIX}_CHH.bw"
else
    echo "Conversion complete. Output files: ${PREFIX}_CG.bw, ${PREFIX}_CHG.bw, ${PREFIX}_CHH.bw"
fi