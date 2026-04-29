/**
 * FollowUpChips — rendered below each response card.
 *
 * Shows 2–3 contextual "what to do next" prompts as clickable pills.
 * Clicking a chip calls onSelect(prompt) which pre-fills and submits the input.
 */
import React from "react";
import { getFollowUpSuggestions } from "../utils/followUpSuggestions";

interface Props {
  tool: string | undefined;
  onSelect: (prompt: string) => void;
}

export const FollowUpChips: React.FC<Props> = ({ tool, onSelect }) => {
  const suggestions = getFollowUpSuggestions(tool);
  if (suggestions.length === 0) return null;

  return (
    <div
      className="mt-3 pt-3"
      style={{ borderTop: "1px solid #F1F5F9" }}
    >
      <div
        className="text-muted mb-2"
        style={{ fontSize: "0.74rem", fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.05em" }}
      >
        What to do next
      </div>
      <div className="d-flex flex-wrap gap-2">
        {suggestions.map((s) => (
          <button
            key={s.prompt}
            onClick={() => onSelect(s.prompt)}
            className="border-0 rounded-pill d-inline-flex align-items-center gap-1"
            style={{
              background: "#F8FAFC",
              border: "1px solid #E2E8F0",
              color: "#334155",
              fontSize: "0.8rem",
              fontWeight: 500,
              padding: "4px 12px",
              cursor: "pointer",
              transition: "background 0.15s, border-color 0.15s",
            }}
            onMouseEnter={(e) => {
              e.currentTarget.style.background = "#EFF6FF";
              e.currentTarget.style.borderColor = "#BFDBFE";
            }}
            onMouseLeave={(e) => {
              e.currentTarget.style.background = "#F8FAFC";
              e.currentTarget.style.borderColor = "#E2E8F0";
            }}
            title={s.prompt}
          >
            <span>{s.icon}</span>
            <span>{s.label}</span>
          </button>
        ))}
      </div>
    </div>
  );
};
