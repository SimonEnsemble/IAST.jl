struct AdsIsoTData
    data::DataFrame
    p_key::String
    l_key::String

    function AdsIsoTData(data::DataFrame, p_key::String, l_key::String)
        # some checks
        @assert p_key in names(data)
        @assert l_key in names(data)

        # sort
        s_data = sort(data, p_key)

        # remove < 0 pressure (unphysical)
        filter!(row -> row[p_key] > 0, s_data)

        return new(s_data, p_key, l_key)
    end
end

# slicing
Base.getindex(d::AdsIsoTData, ids) = AdsIsoTData(d.data[ids, :], d.p_key, d.l_key)
Base.length(d::AdsIsoTData) = nrow(d.data)
Base.lastindex(d::AdsIsoTData) = Base.lastindex(d.data)
Base.size(d::AdsIsoTData) = nrow(d.data)
