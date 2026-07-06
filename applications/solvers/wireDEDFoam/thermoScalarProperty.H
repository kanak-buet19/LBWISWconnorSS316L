#ifndef thermoScalarProperty_H
#define thermoScalarProperty_H

class thermoScalarProperty
{
    const word name_;
    List<Tuple2<scalar, scalar>> table_;
    List<scalar> coefficients_;

    scalar tableValue(const scalar x) const
    {
        if (x <= table_.first().first())
        {
            return table_.first().second();
        }

        if (x >= table_.last().first())
        {
            return table_.last().second();
        }

        for (label i = 1; i < table_.size(); ++i)
        {
            const scalar x0 = table_[i - 1].first();
            const scalar x1 = table_[i].first();

            if (x <= x1)
            {
                const scalar y0 = table_[i - 1].second();
                const scalar y1 = table_[i].second();
                const scalar weight = (x - x0)/(x1 - x0);

                return y0 + weight*(y1 - y0);
            }
        }

        return table_.last().second();
    }

    scalar polynomialValue(const scalar x) const
    {
        scalar result = 0;

        forAllReverse(coefficients_, i)
        {
            result = result*x + coefficients_[i];
        }

        return result;
    }

    void validateTable(const word& tableKeyword) const
    {
        if (table_.empty())
        {
            FatalErrorInFunction
                << "Property table " << tableKeyword << " for " << name_
                << " is empty" << nl
                << exit(FatalError);
        }

        for (label i = 1; i < table_.size(); ++i)
        {
            if (table_[i].first() <= table_[i - 1].first())
            {
                FatalErrorInFunction
                    << "Property table " << tableKeyword << " for " << name_
                    << " must use strictly increasing temperatures. Entry "
                    << i - 1 << " has T=" << table_[i - 1].first()
                    << ", entry " << i << " has T=" << table_[i].first()
                    << nl << exit(FatalError);
            }
        }
    }

public:

    thermoScalarProperty
    (
        const dictionary& dict,
        const word& tableKeyword,
        const word& polynomialKeyword,
        const word& name
    )
    :
        name_(name),
        table_(),
        coefficients_()
    {
        if (dict.found(tableKeyword))
        {
            table_ = List<Tuple2<scalar, scalar>>(dict.lookup(tableKeyword));
            validateTable(tableKeyword);
        }

        if (dict.found(polynomialKeyword))
        {
            coefficients_ = List<scalar>(dict.lookup(polynomialKeyword));
        }

        if (table_.empty() && coefficients_.empty())
        {
            FatalErrorInFunction
                << "No thermophysical property data found for " << name_
                << ". Provide " << tableKeyword << " or "
                << polynomialKeyword << nl
                << exit(FatalError);
        }
    }

    scalar value(const scalar x) const
    {
        if (!table_.empty())
        {
            return tableValue(x);
        }

        return polynomialValue(x);
    }

    bool usesTable() const
    {
        return !table_.empty();
    }
};

#endif
