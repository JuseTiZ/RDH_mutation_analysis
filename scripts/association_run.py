import argparse
import pandas as pd
import statsmodels.api as sm
import statsmodels
import numpy as np
from tqdm import tqdm

def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)

    target_variable = parser.add_argument_group("Target")
    target_variable.add_argument("--covariates", nargs="+", required=True)
    target_variable.add_argument("--var_of_interest", nargs="+", required=True)
    target_variable.add_argument("--mask", nargs="+", default=[])

    phenotype_variable = parser.add_argument_group("Phenotype")
    phenotype_variable.add_argument("-p", "--phenotype", nargs="+", required=True)
    phenotype_variable.add_argument("--sex", nargs="+", default=[])
    phenotype_variable.add_argument("--OLS", action="store_true")

    return parser.parse_args()


def _result_prep_logit(result, var_of_interest_index, debug=False):
    """
    Process result from statsmodels
    :param result: logistic regression result
    :param var_of_interest_index: index of variable of interest
    :return: dataframe with key statistics
    """
    results_as_html = result.summary().tables[0].as_html()
    converged = pd.read_html(results_as_html)[0].iloc[5, 1]
    results_as_html = result.summary().tables[1].as_html()
    res = pd.read_html(results_as_html, header=0, index_col=0)[0]

    if debug:
        print(res)

    p_value = result.pvalues[var_of_interest_index]
    neg_log_p_value = -np.log10(p_value)
    beta = result.params[var_of_interest_index]
    ci_lower = res.iloc[var_of_interest_index]["[0.025"]
    ci_upper = res.iloc[var_of_interest_index]["0.975]"]
    odds_ratio = np.exp(beta)

    return {
        "p_value": p_value,
        "neg_log_p_value": neg_log_p_value,
        "beta": beta,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "odds_ratio": odds_ratio,
        "converged": converged,
    }


def _result_prep_OLS(result, var_of_interest_index, debug=False):
    """
    Process result from statsmodels
    :param result: ordinary least squares result
    :param var_of_interest_index: index of variable of interest
    :return: dataframe with key statistics
    """
    results_as_html = result.summary().tables[0].as_html()
    results_as_html = result.summary().tables[1].as_html()
    res = pd.read_html(results_as_html, header=0, index_col=0)[0]

    if debug:
        print(res)

    p_value = result.pvalues[var_of_interest_index]
    neg_log_p_value = -np.log10(p_value)
    beta = result.params[var_of_interest_index]
    ci_lower = res.iloc[var_of_interest_index]["[0.025"]
    ci_upper = res.iloc[var_of_interest_index]["0.975]"]

    return {
        "p_value": p_value,
        "neg_log_p_value": neg_log_p_value,
        "beta": beta,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
    }


def main():
    # Ignore warnings
    import warnings
    warnings.filterwarnings("ignore")

    args = get_args()
    if not args.mask:
        args.mask = [None] * len(args.var_of_interest)
    else:
        if len(args.mask) != len(args.var_of_interest):
            raise ValueError("Length of mask must be equal to length of var_of_interest")
        args.mask = [None if x == "None" else x for x in args.mask]

    if not args.sex:
        args.sex = [None] * len(args.var_of_interest)
    else:
        if len(args.sex) != len(args.phenotype):
            raise ValueError("Length of sex must be equal to length of phenotype")
        args.sex = [None if x == "None" else x for x in args.sex]

    # Read in dataframe
    df = pd.read_csv(args.input)
    results = []

    # Iterate over variables of interest
    bar = tqdm(zip(args.var_of_interest, args.mask), total=len(args.var_of_interest))
    for variable, mask in bar:
        # Iterate over phenotypes
        for phenotype, sex in zip(args.phenotype, args.sex):

            columns_needed = [variable] + args.covariates + [phenotype]
            if mask is not None:
                columns_needed.append(mask)
            sub_df = df[columns_needed].dropna()
            if mask is not None:
                sub_df = sub_df[sub_df[mask] == 0]

            if sex is not None:
                sub_df = sub_df[sub_df['sex'] == int(sex)]
                using_sex = sex
                sub_df = sub_df.drop(['sex'], axis=1)
                covariates = list(set(args.covariates) - {'sex'})
            else:
                using_sex = "both"
                covariates = args.covariates

            bar.set_description(f"Processing {phenotype} (N {len(sub_df)})")
            
            X = sub_df[covariates + [variable]].to_numpy()
            X = sm.tools.add_constant(X)
            y = sub_df[phenotype].to_numpy()
            # Fit model
            if args.OLS:
                model = sm.OLS(y, X, missing="drop")
            else:
                model = sm.Logit(y, X, missing="drop")

            try:
                result = model.fit(disp=False)
            except (
                np.linalg.linalg.LinAlgError,
                statsmodels.tools.sm_exceptions.PerfectSeparationError,
            ) as err:
                if "Singular matrix" in str(err) or "Perfect separation" in str(err):
                    print(
                        f"Skipping {variable} and {phenotype} due to singular matrix / perfect separation."
                    )
                    continue
                
            # Record results
            base_dict = {
                "variable": variable,
                "phenotype": phenotype,
            }
            stats_dict = (
                _result_prep_OLS(result, -1)
                if args.OLS
                else _result_prep_logit(result, -1)
            )

            if args.OLS:
                result_dict = {**base_dict,
                               "#indivi": len(y),
                               "#indivi_with_variant": sum(sub_df[variable] > 0),
                               **stats_dict,
                               "sex": using_sex,}
            else:
                cases = sum(y == 1)
                controls = sum(y == 0)
                case_no_variant_count = len(
                    sub_df[(sub_df[phenotype] == 1) & (sub_df[variable] == 0)]
                )
                control_no_variant_count = len(
                    sub_df[(sub_df[phenotype] == 0) & (sub_df[variable] == 0)]
                )
                case_have_variant_count = len(
                    sub_df[(sub_df[phenotype] == 1) & (sub_df[variable] > 0)]
                )
                control_have_variant_count = len(
                    sub_df[(sub_df[phenotype] == 0) & (sub_df[variable] > 0)]
                )
                result_dict = {
                    **base_dict,
                    "cases": cases,
                    "controls": controls,
                    **stats_dict,
                    "case_have_variant_count": case_have_variant_count,
                    "control_have_variant_count": control_have_variant_count,
                    "case_no_variant_count": case_no_variant_count,
                    "control_no_variant_count": control_no_variant_count,
                    "sex": using_sex,
                }

            results.append(result_dict)
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()