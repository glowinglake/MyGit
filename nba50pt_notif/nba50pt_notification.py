#!/usr/bin/env python3
"""
NBA 50-Point Performance Notifier

Checks yesterday's NBA games for players who scored 50+ points
and sends macOS desktop notifications if any are found.

Usage:
    python nba50pt_notification.py

Schedule via crontab:
    0 9 * * * /Users/yinyi/tmp/MyGit/nba50pt_notif/nba_venv/bin/python /Users/yinyi/tmp/MyGit/nba50pt_notif/nba50pt_notification.py >> /Users/yinyi/tmp/MyGit/nba50pt_notif/nba50pt.log 2>&1
"""

import subprocess
from datetime import datetime, timedelta
from nba_api.stats.endpoints import ScoreboardV2, BoxScoreTraditionalV3


def send_notification(title: str, message: str) -> None:
    """Send a macOS desktop notification using osascript."""
    # Escape double quotes in the message for AppleScript
    escaped_message = message.replace('"', '\\"')
    escaped_title = title.replace('"', '\\"')
    
    subprocess.run([
        'osascript', '-e',
        f'display notification "{escaped_message}" with title "{escaped_title}"'
    ], check=False)


def get_yesterday_date() -> str:
    """Return yesterday's date in YYYY-MM-DD format."""
    yesterday = datetime.now() - timedelta(days=1)
    return yesterday.strftime('%Y-%m-%d')


def get_game_ids_for_date(game_date: str) -> list[str]:
    """Fetch all game IDs for a given date."""
    print(f"Fetching games for {game_date}...")
    scoreboard = ScoreboardV2(game_date=game_date)
    games_df = scoreboard.get_data_frames()[0]  # GameHeader DataFrame
    
    if games_df.empty:
        print("No games found for this date.")
        return []
    
    game_ids = games_df['GAME_ID'].tolist()
    print(f"Found {len(game_ids)} game(s).")
    return game_ids


def get_50pt_scorers(game_id: str) -> list[dict]:
    """Get all players who scored 50+ points in a given game."""
    scorers = []
    
    try:
        boxscore = BoxScoreTraditionalV3(game_id=game_id)
        player_stats_df = boxscore.player_stats.get_data_frame()
        
        # Filter for players with 50+ points
        high_scorers = player_stats_df[player_stats_df['points'] >= 50]
        
        for _, row in high_scorers.iterrows():
            scorers.append({
                'name': f"{row['firstName']} {row['familyName']}",
                'team': row['teamTricode'],
                'points': int(row['points']),
                'game_id': game_id
            })
    except Exception as e:
        print(f"Error fetching boxscore for game {game_id}: {e}")
    
    return scorers


def main():
    """Main function to check for 50-point performances and notify."""
    yesterday = get_yesterday_date()
    print(f"Checking NBA games from {yesterday} for 50+ point performances...")
    
    # Get all games from yesterday
    game_ids = get_game_ids_for_date(yesterday)
    
    if not game_ids:
        print("No games to check.")
        return
    
    # Check each game for 50+ point scorers
    all_50pt_scorers = []
    for game_id in game_ids:
        scorers = get_50pt_scorers(game_id)
        all_50pt_scorers.extend(scorers)
    
    # Send notifications if any 50+ point performances found
    if all_50pt_scorers:
        print(f"\n🏀 Found {len(all_50pt_scorers)} player(s) with 50+ points!")
        
        for scorer in all_50pt_scorers:
            message = f"{scorer['name']} ({scorer['team']}) scored {scorer['points']} points!"
            print(f"  - {message}")
            send_notification("NBA 50+ Point Alert! 🏀", message)
    else:
        print("\nNo 50+ point performances yesterday.")


if __name__ == '__main__':
    main()

