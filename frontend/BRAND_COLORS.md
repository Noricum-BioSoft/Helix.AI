# 🎨 Helix.AI Brand Colors

This document describes the brand colors extracted from the NORICUM logo and how they're used throughout the frontend.

## Color Palette

### Primary Brand Colors

| Color | Hex Code | Usage | Description |
|-------|----------|-------|-------------|
| **Gold** | `#C8B280` | Primary actions, highlights | Golden-beige from logo banner |
| **Blue** | `#3A60A8` | Secondary actions, links | Blue from "NORICUM" text |
| **White** | `#FFFFFF` | Background | White background |

### Color Variations

#### Gold Variations
- **Gold Light**: `#E8DCC0` - Subtle backgrounds
- **Gold Lighter**: `#F4EDE0` - Very subtle backgrounds
- **Gold Dark**: `#A8966A` - Hover states
- **Gold Darker**: `#8B7A5A` - Active states

#### Blue Variations
- **Blue Light**: `#5A7AB8` - Hover states
- **Blue Lighter**: `#7A9AC8` - Light accents
- **Blue Dark**: `#2A4A88` - Active states
- **Blue Darker**: `#1A3A68` - Deep accents

## Implementation

### Theme Files

1. **`theme.css`** - CSS variables and Bootstrap overrides
   - Defines CSS custom properties for all colors
   - Overrides Bootstrap's default button, badge, and link colors
   - Provides utility classes for brand colors

2. **`theme.ts`** - TypeScript theme object
   - Programmatic access to colors in React components
   - Includes spacing, typography, shadows, and other design tokens
   - Type-safe color references

### Usage Examples

#### CSS Variables
```css
.my-element {
  background-color: var(--brand-gold);
  color: var(--brand-blue);
  border: 1px solid var(--border-gold);
}
```

#### TypeScript/React
```tsx
import { theme } from './theme';

const MyComponent = () => (
  <div style={{ 
    backgroundColor: theme.colors.gold,
    color: theme.colors.blue 
  }}>
    Content
  </div>
);
```

#### Bootstrap Classes
```tsx
<Button variant="primary">  {/* Uses gold */}
<Button variant="secondary"> {/* Uses blue */}
<span className="text-primary">  {/* Gold text */}
<span className="text-secondary"> {/* Blue text */}
```

## Color Application

### Primary (Gold) - `#C8B280`
- **Buttons**: Primary action buttons
- **Backgrounds**: Tips section, empty states
- **Borders**: File upload zone, cards
- **Shadows**: Button shadows, hover effects

### Secondary (Blue) - `#3A60A8`
- **Links**: All clickable links
- **Badges**: Session badges, info badges
- **Text**: Command examples, code snippets
- **Accents**: Hover states, active items

### Backgrounds
- **Gold Subtle** (`#F4EDE0`): Tips panel, upload zone
- **Blue Subtle** (`#E8EDF5`): Hover states, active categories

## Components Updated

All components now use the brand colors:

1. **App.tsx**
   - Tips section (gold background, blue text)
   - Session badge (blue)
   - File upload zone (gold/blue borders)
   - Submit button (gold with gold shadow)

2. **WelcomeScreen.tsx**
   - Primary button (gold)
   - Example buttons (gold outline)

3. **EmptyState.tsx**
   - Action cards (gold borders on hover)
   - Upload button (gold outline)
   - Example button (blue outline)

4. **ExampleCommandsPanel.tsx**
   - Header (gold background)
   - Categories (gold/blue borders)
   - Command text (blue)
   - Hover states (blue backgrounds)

## Bootstrap Overrides

The following Bootstrap classes are overridden:

- `.btn-primary` → Uses gold
- `.btn-secondary` → Uses blue
- `.btn-outline-primary` → Gold outline
- `.btn-outline-secondary` → Blue outline
- `.text-primary` → Gold text
- `.text-secondary` → Blue text
- `.bg-primary` → Gold background
- `.bg-secondary` → Blue background
- `.badge.bg-info` → Blue badge
- `a` (links) → Blue links

## Accessibility

- **Gold on White**: Good contrast for primary actions
- **Blue on White**: Excellent contrast for links and text
- **Gold Hover**: Darker gold for better hover feedback
- **Blue Hover**: Darker blue for better hover feedback

## Color Psychology

- **Gold**: Represents value, quality, and scientific excellence
- **Blue**: Represents trust, reliability, and technology
- **Combination**: Creates a professional, scientific, and trustworthy appearance

## Future Enhancements

- Dark mode support with adjusted brand colors
- Additional color variations for different contexts
- Gradient support using brand colors
- Color palette expansion for data visualization

---

**Last Updated:** $(date)  
**Source:** NORICUM Logo  
**Implementation:** CSS Variables + TypeScript Theme Object



