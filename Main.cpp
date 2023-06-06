void solve_distribution_problem(vector<double> b,vector<vector<vector<double>>>
red,→ production_tables,vector<vector<double>> alphas, vector<vector<double>>
red,→ subs_costs, vector<vector<double>> main_costs)
{
size_t subs = production_tables.size();
size_t supply_size = b.size();
vector<map<vector<double>, double>> subs_productions;
auto possible_supplies = generate_possible_supplies(b);
for (size_t sub = 0; sub < subs; ++sub)
{
for (auto& possible_supply : possible_supplies)
{
subs_productions[sub][possible_supply] = INT_MIN;
}
}
for (auto& possible_supply : possible_supplies)
{
subs_productions[0][possible_supply] =
scalar_product(solve_linear_programming(production_tables[0], sum_vectors(
red,→ alphas[0], possible_supply), subs_costs[0]), main_costs[0]);
}
for (size_t i = 1; i < subs; ++i)
{
for (auto possible_supply : possible_supplies)
{
auto possible_supplies_for_current_sub = generate_possible_supplies(
red,→ possible_supply);
double production = 0;
for (auto possible_supply_for_current_sub : possible_supplies_for_current_sub)
{
auto left_supplies = subtract_vectors(possible_supply,
red,→ possible_supply_for_current_sub);
production = max(scalar_product(solve_linear_programming(
red,→ production_tables[i], sum_vectors(alphas[i],
red,→ possible_supply_for_current_sub), subs_costs[i]), main_costs[i]) +
subs_productions[i − 1][left_supplies], production);
subs_productions[i][possible_supply] = production;
}
}
}

int main()
{
     solve_distribution_problem(random_matrix);
}
