/**
 * Verbose UI logging is off in production unless explicitly enabled.
 *
 * Enable detailed logs:
 * - Append `?debug=1` or `?helix_debug=1` to the URL (then reload if needed)
 * - Run `localStorage.setItem('HELIX_DEBUG', '1')` in the devtools console and reload
 * - Build with `VITE_HELIX_DEBUG=true`
 * - Vite dev server: logging is on when `import.meta.env.DEV` is true
 */

function readDebugFlag(): boolean {
  if (import.meta.env.DEV) return true;
  if (import.meta.env.VITE_HELIX_DEBUG === 'true') return true;
  if (typeof window === 'undefined') return false;
  try {
    const q = new URLSearchParams(window.location.search);
    if (q.get('debug') === '1' || q.get('helix_debug') === '1') return true;
    if (window.localStorage?.getItem('HELIX_DEBUG') === '1') return true;
  } catch {
    /* ignore */
  }
  return false;
}

const DEBUG = readDebugFlag();

export function isHelixDebugEnabled(): boolean {
  return DEBUG;
}

/** Verbose trace; no-op in production unless debug is enabled (see module doc). */
export function debugLog(...args: unknown[]): void {
  if (DEBUG) console.log('[Helix]', ...args);
}

export function debugWarn(...args: unknown[]): void {
  if (DEBUG) console.warn('[Helix]', ...args);
}
