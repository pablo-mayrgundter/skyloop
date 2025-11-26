module.exports = {
  root: true,
  env: {
    browser: true,
    es2021: true,
    node: true
  },
  parser: '@typescript-eslint/parser',
  plugins: ['@typescript-eslint', 'react'],
  extends: [
    'eslint:recommended',
    'plugin:react/recommended',
    'google'
  ],
  parserOptions: {
    ecmaVersion: 'latest',
    sourceType: 'module'
  },
  settings: {
    react: {
      version: 'detect'
    }
  },
  rules: {
    'semi': ['error', 'never'],
    'quote-props': ['error', 'as-needed'],
    'indent': ['error', 2],
    'max-len': ['error', 120],
    'object-curly-spacing': ['error', 'always'],
    'react/react-in-jsx-scope': 'off'
  }
}
