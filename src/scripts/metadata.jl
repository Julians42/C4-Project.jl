# scraping and adding instrument response for SCEDC_data
using AWS, AWSS3, SeisIO, SeisNoise, Dates


function yyyyjjj2date(yearday::String)
    @assert occursin(r"[1-2][0-9][0-9][0-9][0-3][0-9][0-9]",yearday)
    yint = parse(Int,yearday[1:4])
    dint = parse(Int,yearday[5:end])
    @assert dint <= 366 "Input day must be less than or equal to 366"
    return DateTime(yint) + Day(dint-1)
end

function read_resp(file::String,XMLDIR::String)
    s = yyyyjjj2date(file[end-9:end-3])
	t = s + Day(1)
	s = Dates.format(s, "yyyy-mm-dd HH:MM:SS")
	t = Dates.format(t, "yyyy-mm-dd HH:MM:SS")
	net = basename(file)[1:2]
	sta = split(basename(file),"_")[1][3:end]
	instpath = joinpath(XMLDIR,net * '_' * sta * ".xml" )
    return read_meta("sxml",instpath,s=s,t=t)
end

function XML_download(aws,XMLDIR)
    if !isdir(XMLDIR)
        mkpath(XMLDIR)
    end
    req = collect(s3_list_objects(aws,"scedc-pds","FDSNstationXML/CI/"))
    xmlin = [r["Key"] for r in req]
    xmlout = joinpath.(XMLDIR,basename.(xmlin))
    for ii = 1:length(xmlin)
        s3_get_file(aws,"scedc-pds",xmlin[ii],xmlout[ii])
    end
    return nothing
end
XML_download(XMLDIR) = XML_download(global_aws_config(region="us-west-2"),XMLDIR)



XMLDIR = "/Volumes/T7/seis_data/XML/"
aws = AWS.AWSConfig(region = "us-west-2")

# download metadata 
XML_download(aws, XMLDIR)

# sample seisdata
S = read_data("mseed", "/Users/julianschmitt/Downloads/CIADO__BHE___2018001.ms")
function add_response(S::SeisData, file::String)
    resp = read_resp(file, XMLDIR) 
    found_responses = resp[findfirst(x -> x == S.id[1], resp.id)] # filter response
    S.resp[1] = found_responses.resp[1] # add response
    S.loc[1] = found_responses.loc[1]
    S.gain[1] = found_responses.gain[1]
end

R = remove_resp(S)
