// Helix.AI Brand Colors - Extracted from NORICUM Logo
export const brandColors = {
  // Primary Brand Colors
  gold: '#C8B280',           // Golden-beige from logo banner
  blue: '#3A60A8',           // Blue from "NORICUM" text
  white: '#FFFFFF',          // White background
  
  // Gold Variations
  goldLight: '#E8DCC0',
  goldLighter: '#F4EDE0',
  goldDark: '#A8966A',
  goldDarker: '#8B7A5A',
  
  // Blue Variations
  blueLight: '#5A7AB8',
  blueLighter: '#7A9AC8',
  blueDark: '#2A4A88',
  blueDarker: '#1A3A68',
  
  // Semantic Colors
  primary: '#C8B280',
  secondary: '#3A60A8',
  
  // Text Colors
  textPrimary: '#2C3E50',
  textSecondary: '#6C757D',
  textOnGold: '#2C3E50',
  textOnBlue: '#FFFFFF',
  
  // Background Colors
  bgPrimary: '#FFFFFF',
  bgSecondary: '#F8F9FA',
  bgGoldSubtle: '#F4EDE0',
  bgBlueSubtle: '#E8EDF5',
  
  // Border Colors
  borderColor: '#DEE2E6',
  borderGold: '#C8B280',
  borderBlue: '#3A60A8',
};

export const theme = {
  colors: brandColors,
  
  // Spacing
  spacing: {
    xs: '0.25rem',
    sm: '0.5rem',
    md: '1rem',
    lg: '1.5rem',
    xl: '2rem',
    xxl: '3rem',
  },
  
  // Typography
  typography: {
    fontFamily: {
      sans: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
      mono: 'Consolas, "Liberation Mono", Menlo, Courier, monospace',
    },
    fontSize: {
      xs: '0.75rem',
      sm: '0.875rem',
      base: '1rem',
      lg: '1.125rem',
      xl: '1.25rem',
      '2xl': '1.5rem',
      '3xl': '1.875rem',
      '4xl': '2.25rem',
    },
  },
  
  // Border Radius
  borderRadius: {
    sm: '0.125rem',
    md: '0.25rem',
    lg: '0.5rem',
    xl: '0.75rem',
    full: '9999px',
  },
  
  // Shadows
  shadows: {
    sm: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',
    md: '0 4px 6px -1px rgba(0, 0, 0, 0.1)',
    lg: '0 10px 15px -3px rgba(0, 0, 0, 0.1)',
    gold: '0 4px 8px rgba(200, 178, 128, 0.3)',
    blue: '0 4px 8px rgba(58, 96, 168, 0.3)',
  },
};

export default theme;



