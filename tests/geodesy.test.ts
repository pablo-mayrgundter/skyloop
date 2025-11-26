import { EARTH_RENDER_RADIUS, LOS_ANGELES, SAN_FRANCISCO, greatCirclePoint, interpolateGreatCircle } from '../src/config/simulation'

describe('great-circle helpers', () => {
  it('returns normalized radius', () => {
    const la = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon, EARTH_RENDER_RADIUS)
    expect(Math.abs(la.length() - EARTH_RENDER_RADIUS)).toBeLessThan(1e-6)
  })

  it('interpolates between endpoints', () => {
    const mid = interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, 0.5)
    expect(mid.length()).toBeCloseTo(EARTH_RENDER_RADIUS)
    const start = greatCirclePoint(LOS_ANGELES.lat, LOS_ANGELES.lon)
    const end = greatCirclePoint(SAN_FRANCISCO.lat, SAN_FRANCISCO.lon)
    expect(mid.distanceTo(start)).toBeGreaterThan(0)
    expect(mid.distanceTo(end)).toBeGreaterThan(0)
  })
})
