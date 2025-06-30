import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64

from smolagents import CodeAgent, tool


class DataAnalysisAgent(CodeAgent):
    """Extended CodeAgent with dataset awareness"""

    def __init__(self, dataset: pd.DataFrame, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._dataset = dataset

    @property
    def dataset(self) -> pd.DataFrame:
        """Access the stored dataset"""
        return self._dataset

    def run(self, prompt: str) -> str:
        """Override run method to include dataset context"""
        dataset_info = f"""
        Dataset Shape: {self.dataset.shape}
        Columns: {', '.join(self.dataset.columns)}
        Data Types: {self.dataset.dtypes.to_dict()}
        """
        enhanced_prompt = f"""
        Analyze the following dataset:
        {dataset_info}

        Task: {prompt}

        Use the provided tools to analyze this specific dataset and return detailed results.
        """

        print(f'[DataAnalysisAgent.run]: {enhanced_prompt}')
        return super().run(enhanced_prompt)


@tool
def analyze_basic_stats(data: pd.DataFrame) -> str:
    """Calculate basic statistical measures for numerical columns in the dataset.

    This function computes fundamental statistical metrics including mean, median,
    standard deviation, skewness, and counts of missing values for all numerical
    columns in the provided DataFrame.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least one numerical column for meaningful analysis.

    Returns:
        str: A string containing formatted basic statistics for each numerical column,
            including mean, median, standard deviation, skewness, and missing value counts.
    """
    # Access dataset from agent if no data provided
    if data is None:
        data = tool.agent.dataset

    stats = {}
    numeric_cols = data.select_dtypes(include=[np.number]).columns

    for col in numeric_cols:
        stats[col] = {
            'mean': float(data[col].mean()),
            'median': float(data[col].median()),
            'std': float(data[col].std()),
            'skew': float(data[col].skew()),
            'missing': int(data[col].isnull().sum())
        }

    return str(stats)


@tool
def generate_correlation_matrix(data: pd.DataFrame) -> str:
    """Generate a visual correlation matrix for numerical columns in the dataset.

    This function creates a heatmap visualization showing the correlations between
    all numerical columns in the dataset. The correlation values are displayed
    using a color-coded matrix for easy interpretation.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least two numerical columns for correlation analysis.

    Returns:
        str: A base64 encoded string representing the correlation matrix plot image,
            which can be displayed in a web interface or saved as an image file.
    """
    # Access dataset from agent if no data provided
    if data is None:
        data = tool.agent.dataset

    numeric_data = data.select_dtypes(include=[np.number])

    plt.figure(figsize=(10, 8))
    sns.heatmap(numeric_data.corr(), annot=True, cmap='coolwarm')
    plt.title('Correlation Matrix')

    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    return base64.b64encode(buf.getvalue()).decode()


@tool
def analyze_categorical_columns(data: pd.DataFrame) -> str:
    """Analyze categorical columns in the dataset for distribution and frequencies.

    This function examines categorical columns to identify unique values, top categories,
    and missing value counts, providing insights into the categorical data distribution.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least one categorical column for meaningful analysis.

    Returns:
        str: A string containing formatted analysis results for each categorical column,
            including unique value counts, top categories, and missing value counts.
    """
    # Access dataset from agent if no data provided
    if data is None:
        data = tool.agent.dataset

    categorical_cols = data.select_dtypes(include=['object', 'category']).columns
    analysis = {}

    for col in categorical_cols:
        analysis[col] = {
            'unique_values': int(data[col].nunique()),
            'top_categories': data[col].value_counts().head(5).to_dict(),
            'missing': int(data[col].isnull().sum())
        }

    return str(analysis)


@tool
def suggest_features(data: pd.DataFrame) -> str:
    """Suggest potential feature engineering steps based on data characteristics.

    This function analyzes the dataset's structure and statistical properties to
    recommend possible feature engineering steps that could improve model performance.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            can contain both numerical and categorical columns.

    Returns:
        str: A string containing suggestions for feature engineering based on
            the characteristics of the input data.
    """
    # Access dataset from agent if no data provided
    if data is None:
        data = tool.agent.dataset

    suggestions = []
    numeric_cols = data.select_dtypes(include=[np.number]).columns
    categorical_cols = data.select_dtypes(include=['object', 'category']).columns

    if len(numeric_cols) >= 2:
        suggestions.append("Consider creating interaction terms between numerical features")

    if len(categorical_cols) > 0:
        suggestions.append("Consider one-hot encoding for categorical variables")

    for col in numeric_cols:
        if data[col].skew() > 1 or data[col].skew() < -1:
            suggestions.append(f"Consider log transformation for {col} due to skewness")

    return '\n'.join(suggestions)