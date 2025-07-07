import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io
import base64
from typing import Dict, Any, List

# from smolagents import CodeAgent, tool


class DataAnalysisAgent:
    """Data analysis agent for bioinformatics datasets"""

    def __init__(self, dataset: pd.DataFrame):
        self._dataset = dataset

    @property
    def dataset(self) -> pd.DataFrame:
        """Access the stored dataset"""
        return self._dataset

    def run(self, prompt: str) -> str:
        """Run analysis on the dataset"""
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
        return enhanced_prompt


def analyze_basic_stats(data: pd.DataFrame) -> Dict[str, Any]:
    """Calculate basic statistical measures for numerical columns in the dataset.

    This function computes fundamental statistical metrics including mean, median,
    standard deviation, skewness, and counts of missing values for all numerical
    columns in the provided DataFrame.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least one numerical column for meaningful analysis.

    Returns:
        Dict: A dictionary containing formatted basic statistics for each numerical column,
            including mean, median, standard deviation, skewness, and missing value counts.
    """
    if data is None or data.empty:
        return {
            "text": "Error: No data provided or data is empty",
            "statistics": {},
            "plot": {"data": [], "layout": {"title": "No Data"}}
        }

    try:
        stats = {}
        numeric_cols = data.select_dtypes(include=[np.number]).columns

        if len(numeric_cols) == 0:
            return {
                "text": "No numerical columns found in the dataset",
                "statistics": {},
                "plot": {"data": [], "layout": {"title": "No Numerical Data"}}
            }

        for col in numeric_cols:
            col_data = data[col].dropna()
            if len(col_data) > 0:
                stats[col] = {
                    'count': int(len(col_data)),
                    'mean': float(col_data.mean()),
                    'median': float(col_data.median()),
                    'std': float(col_data.std()),
                    'min': float(col_data.min()),
                    'max': float(col_data.max()),
                    'skew': float(col_data.skew()),
                    'missing': int(data[col].isnull().sum())
                }
            else:
                stats[col] = {
                    'count': 0,
                    'mean': 0.0,
                    'median': 0.0,
                    'std': 0.0,
                    'min': 0.0,
                    'max': 0.0,
                    'skew': 0.0,
                    'missing': int(data[col].isnull().sum())
                }

        # Create visualization
        if len(numeric_cols) > 0:
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle('Basic Statistics Overview')

            # Histogram of first numerical column
            if len(numeric_cols) > 0:
                col = numeric_cols[0]
                axes[0, 0].hist(data[col].dropna(), bins=20, alpha=0.7)
                axes[0, 0].set_title(f'Distribution of {col}')
                axes[0, 0].set_xlabel(col)
                axes[0, 0].set_ylabel('Frequency')

            # Box plot
            if len(numeric_cols) > 0:
                data[numeric_cols].boxplot(ax=axes[0, 1])
                axes[0, 1].set_title('Box Plot of Numerical Columns')
                axes[0, 1].tick_params(axis='x', rotation=45)

            # Missing values
            missing_counts = [stats[col]['missing'] for col in numeric_cols]
            axes[1, 0].bar(numeric_cols, missing_counts)
            axes[1, 0].set_title('Missing Values by Column')
            axes[1, 0].set_ylabel('Missing Count')
            axes[1, 0].tick_params(axis='x', rotation=45)

            # Summary statistics
            summary_data = []
            for col in numeric_cols:
                summary_data.append([
                    stats[col]['mean'],
                    stats[col]['median'],
                    stats[col]['std']
                ])
            
            axes[1, 1].table(cellText=summary_data,
                            rowLabels=numeric_cols,
                            colLabels=['Mean', 'Median', 'Std'],
                            cellLoc='center',
                            loc='center')
            axes[1, 1].set_title('Summary Statistics')
            axes[1, 1].axis('off')

            plt.tight_layout()
            
            # Convert plot to base64
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            plt.close()
            plot_data = base64.b64encode(buf.getvalue()).decode()
        else:
            plot_data = ""

        # Create result text
        result_text = f"Basic Statistics Analysis\n\n"
        result_text += f"Dataset shape: {data.shape}\n"
        result_text += f"Numerical columns: {len(numeric_cols)}\n\n"
        
        for col in numeric_cols:
            result_text += f"{col}:\n"
            result_text += f"  Count: {stats[col]['count']}\n"
            result_text += f"  Mean: {stats[col]['mean']:.4f}\n"
            result_text += f"  Median: {stats[col]['median']:.4f}\n"
            result_text += f"  Std: {stats[col]['std']:.4f}\n"
            result_text += f"  Min: {stats[col]['min']:.4f}\n"
            result_text += f"  Max: {stats[col]['max']:.4f}\n"
            result_text += f"  Skew: {stats[col]['skew']:.4f}\n"
            result_text += f"  Missing: {stats[col]['missing']}\n\n"

        return {
            "text": result_text,
            "statistics": stats,
            "plot": {
                "data": [{"type": "image", "source": f"data:image/png;base64,{plot_data}"}],
                "layout": {"title": "Basic Statistics Overview"}
            }
        }

    except Exception as e:
        return {
            "text": f"Error in basic statistics analysis: {str(e)}",
            "statistics": {},
            "plot": {"data": [], "layout": {"title": "Error"}}
        }


def generate_correlation_matrix(data: pd.DataFrame) -> Dict[str, Any]:
    """Generate a visual correlation matrix for numerical columns in the dataset.

    This function creates a heatmap visualization showing the correlations between
    all numerical columns in the dataset. The correlation values are displayed
    using a color-coded matrix for easy interpretation.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least two numerical columns for correlation analysis.

    Returns:
        Dict: A dictionary containing the correlation matrix plot and statistics.
    """
    if data is None or data.empty:
        return {
            "text": "Error: No data provided or data is empty",
            "correlation_matrix": {},
            "plot": {"data": [], "layout": {"title": "No Data"}}
        }

    try:
        numeric_data = data.select_dtypes(include=[np.number])

        if len(numeric_data.columns) < 2:
            return {
                "text": "Error: At least 2 numerical columns are required for correlation analysis",
                "correlation_matrix": {},
                "plot": {"data": [], "layout": {"title": "Insufficient Data"}}
            }

        # Calculate correlation matrix
        corr_matrix = numeric_data.corr()

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0,
                   square=True, linewidths=0.5)
        plt.title('Correlation Matrix')
        plt.tight_layout()

        # Convert plot to base64
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        plt.close()
        plot_data = base64.b64encode(buf.getvalue()).decode()

        # Create result text
        result_text = f"Correlation Matrix Analysis\n\n"
        result_text += f"Number of numerical columns: {len(numeric_data.columns)}\n\n"
        
        # Find strongest correlations
        corr_pairs = []
        for i in range(len(corr_matrix.columns)):
            for j in range(i+1, len(corr_matrix.columns)):
                col1 = corr_matrix.columns[i]
                col2 = corr_matrix.columns[j]
                corr_value = corr_matrix.iloc[i, j]
                corr_pairs.append((col1, col2, corr_value))

        # Sort by absolute correlation value
        corr_pairs.sort(key=lambda x: abs(x[2]), reverse=True)
        
        result_text += "Top correlations:\n"
        for col1, col2, corr_value in corr_pairs[:5]:
            result_text += f"  {col1} vs {col2}: {corr_value:.4f}\n"

        return {
            "text": result_text,
            "correlation_matrix": corr_matrix.to_dict(),
            "plot": {
                "data": [{"type": "image", "source": f"data:image/png;base64,{plot_data}"}],
                "layout": {"title": "Correlation Matrix"}
            }
        }

    except Exception as e:
        return {
            "text": f"Error in correlation analysis: {str(e)}",
            "correlation_matrix": {},
            "plot": {"data": [], "layout": {"title": "Error"}}
        }


def analyze_categorical_columns(data: pd.DataFrame) -> Dict[str, Any]:
    """Analyze categorical columns in the dataset for distribution and frequencies.

    This function examines categorical columns to identify unique values, top categories,
    and missing value counts, providing insights into the categorical data distribution.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            should contain at least one categorical column for meaningful analysis.

    Returns:
        Dict: A dictionary containing formatted analysis results for each categorical column,
            including unique value counts, top categories, and missing value counts.
    """
    if data is None or data.empty:
        return {
            "text": "Error: No data provided or data is empty",
            "categorical_analysis": {},
            "plot": {"data": [], "layout": {"title": "No Data"}}
        }

    try:
        categorical_cols = data.select_dtypes(include=['object', 'category']).columns
        analysis = {}

        if len(categorical_cols) == 0:
            return {
                "text": "No categorical columns found in the dataset",
                "categorical_analysis": {},
                "plot": {"data": [], "layout": {"title": "No Categorical Data"}}
            }

        for col in categorical_cols:
            analysis[col] = {
                'unique_values': int(data[col].nunique()),
                'top_categories': data[col].value_counts().head(10).to_dict(),
                'missing': int(data[col].isnull().sum()),
                'total_count': int(len(data[col]))
            }

        # Create visualization
        if len(categorical_cols) > 0:
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle('Categorical Data Analysis')

            # Plot 1: Unique values per column
            unique_counts = [analysis[col]['unique_values'] for col in categorical_cols]
            axes[0, 0].bar(categorical_cols, unique_counts)
            axes[0, 0].set_title('Unique Values per Column')
            axes[0, 0].set_ylabel('Unique Count')
            axes[0, 0].tick_params(axis='x', rotation=45)

            # Plot 2: Missing values
            missing_counts = [analysis[col]['missing'] for col in categorical_cols]
            axes[0, 1].bar(categorical_cols, missing_counts)
            axes[0, 1].set_title('Missing Values per Column')
            axes[0, 1].set_ylabel('Missing Count')
            axes[0, 1].tick_params(axis='x', rotation=45)

            # Plot 3: Top categories for first categorical column
            if len(categorical_cols) > 0:
                col = categorical_cols[0]
                top_cats = data[col].value_counts().head(10)
                axes[1, 0].bar(range(len(top_cats)), top_cats.values)
                axes[1, 0].set_title(f'Top Categories in {col}')
                axes[1, 0].set_ylabel('Count')
                axes[1, 0].set_xticks(range(len(top_cats)))
                axes[1, 0].set_xticklabels(top_cats.index, rotation=45)

            # Plot 4: Summary table
            summary_data = []
            for col in categorical_cols:
                summary_data.append([
                    analysis[col]['unique_values'],
                    analysis[col]['missing'],
                    analysis[col]['total_count']
                ])
            
            axes[1, 1].table(cellText=summary_data,
                            rowLabels=categorical_cols,
                            colLabels=['Unique', 'Missing', 'Total'],
                            cellLoc='center',
                            loc='center')
            axes[1, 1].set_title('Categorical Summary')
            axes[1, 1].axis('off')

            plt.tight_layout()
            
            # Convert plot to base64
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            plt.close()
            plot_data = base64.b64encode(buf.getvalue()).decode()
        else:
            plot_data = ""

        # Create result text
        result_text = f"Categorical Data Analysis\n\n"
        result_text += f"Number of categorical columns: {len(categorical_cols)}\n\n"
        
        for col in categorical_cols:
            result_text += f"{col}:\n"
            result_text += f"  Unique values: {analysis[col]['unique_values']}\n"
            result_text += f"  Missing values: {analysis[col]['missing']}\n"
            result_text += f"  Total count: {analysis[col]['total_count']}\n"
            result_text += f"  Top categories:\n"
            for category, count in list(analysis[col]['top_categories'].items())[:5]:
                result_text += f"    {category}: {count}\n"
            result_text += "\n"

        return {
            "text": result_text,
            "categorical_analysis": analysis,
            "plot": {
                "data": [{"type": "image", "source": f"data:image/png;base64,{plot_data}"}],
                "layout": {"title": "Categorical Data Analysis"}
            }
        }

    except Exception as e:
        return {
            "text": f"Error in categorical analysis: {str(e)}",
            "categorical_analysis": {},
            "plot": {"data": [], "layout": {"title": "Error"}}
        }


def suggest_features(data: pd.DataFrame) -> Dict[str, Any]:
    """Suggest potential feature engineering steps based on data characteristics.

    This function analyzes the dataset's structure and statistical properties to
    recommend possible feature engineering steps that could improve model performance.

    Args:
        data: A pandas DataFrame containing the dataset to analyze. The DataFrame
            can contain both numerical and categorical columns.

    Returns:
        Dict: A dictionary containing suggestions for feature engineering based on
            the characteristics of the input data.
    """
    if data is None or data.empty:
        return {
            "text": "Error: No data provided or data is empty",
            "suggestions": [],
            "plot": {"data": [], "layout": {"title": "No Data"}}
        }

    try:
        suggestions = []
        numeric_cols = data.select_dtypes(include=[np.number]).columns
        categorical_cols = data.select_dtypes(include=['object', 'category']).columns

        # Analyze numerical features
        for col in numeric_cols:
            col_data = data[col].dropna()
            if len(col_data) > 0:
                skewness = col_data.skew()
                if abs(skewness) > 1:
                    suggestions.append(f"Consider log transformation for '{col}' (skewness: {skewness:.2f})")
                
                if col_data.min() < 0:
                    suggestions.append(f"Consider absolute value transformation for '{col}' (contains negative values)")
                
                if col_data.std() > col_data.mean():
                    suggestions.append(f"Consider normalization for '{col}' (high variance)")

        # Analyze categorical features
        for col in categorical_cols:
            unique_count = data[col].nunique()
            if unique_count > 10:
                suggestions.append(f"Consider grouping rare categories in '{col}' ({unique_count} unique values)")
            elif unique_count == 2:
                suggestions.append(f"Consider binary encoding for '{col}' (binary categorical)")

        # Interaction suggestions
        if len(numeric_cols) >= 2:
            suggestions.append("Consider creating interaction terms between numerical features")
        
        if len(categorical_cols) > 0 and len(numeric_cols) > 0:
            suggestions.append("Consider creating interaction terms between categorical and numerical features")

        # Missing value suggestions
        missing_cols = data.columns[data.isnull().sum() > 0]
        if len(missing_cols) > 0:
            suggestions.append(f"Consider imputation strategies for columns with missing values: {list(missing_cols)}")

        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Feature Engineering Suggestions')

        # Plot 1: Data types distribution
        dtype_counts = data.dtypes.value_counts()
        axes[0, 0].pie(dtype_counts.values, labels=dtype_counts.index, autopct='%1.1f%%')
        axes[0, 0].set_title('Data Types Distribution')

        # Plot 2: Missing values
        missing_counts = data.isnull().sum()
        missing_counts = missing_counts[missing_counts > 0]
        if len(missing_counts) > 0:
            axes[0, 1].bar(missing_counts.index, missing_counts.values)
            axes[0, 1].set_title('Missing Values')
            axes[0, 1].tick_params(axis='x', rotation=45)
        else:
            axes[0, 1].text(0.5, 0.5, 'No Missing Values', ha='center', va='center')
            axes[0, 1].set_title('Missing Values')

        # Plot 3: Numerical features skewness
        if len(numeric_cols) > 0:
            skewness_values = [data[col].skew() for col in numeric_cols]
            axes[1, 0].bar(numeric_cols, skewness_values)
            axes[1, 0].set_title('Skewness of Numerical Features')
            axes[1, 0].tick_params(axis='x', rotation=45)
            axes[1, 0].axhline(y=1, color='r', linestyle='--', alpha=0.7)
            axes[1, 0].axhline(y=-1, color='r', linestyle='--', alpha=0.7)

        # Plot 4: Categorical features cardinality
        if len(categorical_cols) > 0:
            cardinality = [data[col].nunique() for col in categorical_cols]
            axes[1, 1].bar(categorical_cols, cardinality)
            axes[1, 1].set_title('Cardinality of Categorical Features')
            axes[1, 1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        
        # Convert plot to base64
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
        plt.close()
        plot_data = base64.b64encode(buf.getvalue()).decode()

        # Create result text
        result_text = f"Feature Engineering Suggestions\n\n"
        result_text += f"Dataset shape: {data.shape}\n"
        result_text += f"Numerical columns: {len(numeric_cols)}\n"
        result_text += f"Categorical columns: {len(categorical_cols)}\n\n"
        
        if suggestions:
            result_text += "Suggestions:\n"
            for i, suggestion in enumerate(suggestions, 1):
                result_text += f"{i}. {suggestion}\n"
        else:
            result_text += "No specific feature engineering suggestions at this time.\n"

        return {
            "text": result_text,
            "suggestions": suggestions,
            "plot": {
                "data": [{"type": "image", "source": f"data:image/png;base64,{plot_data}"}],
                "layout": {"title": "Feature Engineering Analysis"}
            }
        }

    except Exception as e:
        return {
            "text": f"Error in feature engineering analysis: {str(e)}",
            "suggestions": [],
            "plot": {"data": [], "layout": {"title": "Error"}}
        }